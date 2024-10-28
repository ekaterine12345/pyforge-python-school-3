from airflow import DAG
from airflow.operators.python import PythonOperator
from datetime import datetime, timedelta
import pandas as pd
from rdkit.Chem import AllChem, Descriptors
from io import BytesIO
import asyncio
from airflow.operators.python import PythonOperator
from airflow.providers.postgres.hooks.postgres import PostgresHook
from airflow.providers.amazon.aws.hooks.s3 import S3Hook
import openpyxl


bucket_name = "airflow"

# Default arguments for the DAG
default_args = {
    'owner': 'coder2j',
    'retries': 5,
    'retry_delay': timedelta(minutes=10)
}


def extract_data():
    hook = PostgresHook(postgres_conn_id="postgres_localhost")
    conn = hook.get_conn()
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM molecules")
    
    # Then fetch the results
    molecules = cursor.fetchall()

    today = datetime.today().date()
    print(molecules)
    # molecules = asyncio.run(MoleculeDao.get_all_mols())
    molecules_dicts = [
        {"identifier": mol[0], "name": mol[1], "smiles": mol[2], "created_at": mol[3]}
        for mol in molecules
    ]

    return molecules_dicts


def transform_data(ti):
    molecules = ti.xcom_pull(task_ids='extract_data')

    mols_df = pd.DataFrame(molecules)

        # check id molecule is okay
    def is_valid_smiles(smiles):
        try:
            mol = AllChem.MolFromSmiles(smiles)
            return mol is not None
        except ValueError:
            return False

    mols_df['is_valid_smiles'] = mols_df['smiles'].apply(is_valid_smiles)
    mols_df = mols_df.loc[mols_df['is_valid_smiles']].copy()

    mols_df['mol'] = mols_df['smiles'].apply(AllChem.MolFromSmiles)

    mol_props_funcs = {
            'Molecular weight': lambda mol: Descriptors.MolWt(mol),
            'logP': lambda mol: Descriptors.MolLogP(mol),
            'H Acceptors': lambda mol: Descriptors.NumHAcceptors(mol),
            'H Donors': lambda mol: Descriptors.NumHDonors(mol),
            'TPSA': lambda mol: Descriptors.TPSA(mol)
    }

    mol_props_to_compute = list(mol_props_funcs.keys())

    mols_df[mol_props_to_compute] = mols_df.apply(
        lambda row: [mol_props_funcs[prop](row['mol']) for prop in mol_props_to_compute],
        axis=1,
        result_type='expand'
    )

    mols_df['Lipinski pass'] = (
        (mols_df['Molecular weight'] < 500) &
        (mols_df['logP'] < 5) &
        (mols_df['H Acceptors'] < 10) &
        (mols_df['H Donors'] < 5)
    )

    mols_df.drop(columns=['mol', 'is_valid_smiles'], inplace=True)
    return mols_df

def save_data(ti):
    mols_df = ti.xcom_pull(task_ids='transform_data')

        # Save to my Excel in memory
    output = BytesIO()
    mols_df.to_excel(output, index=False)
    output.seek(0)

    # to upload the file to MinIO
    file_name = f"molecule_properties_{datetime.now().strftime('%Y-%m-%d')}.xlsx"

    s3_hook = S3Hook(aws_conn_id="minio_conn") 

    # Upload to MinIO bucket
    s3_hook.load_bytes(
        output.getvalue(),
        key=file_name,
        bucket_name="airflow", 
        replace=True
    )

    print(f"File {file_name} successfully uploaded to MinIO bucket 'airflow'")


# Define the DAG
with DAG(
    'molecule_processing_dag_03', 
    default_args=default_args,
    start_date=datetime(2024, 10, 28),
    schedule_interval='@daily'
) as dag:

    # Task 1: Extract Data
    task1 = PythonOperator(
        task_id="extract_data",
        python_callable=extract_data
    )
    

    # Task 2: Transform Data
    task2 = PythonOperator(
        task_id="transform_data",
        python_callable=transform_data
    )

    # Task 3: Save Data and Upload to MinIO
    task3 = PythonOperator(
        task_id="save_data",
        python_callable=save_data
    )

    task1 >> task2 >> task3
