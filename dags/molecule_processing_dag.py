from airflow import DAG
from airflow.operators.python import PythonOperator
from datetime import datetime
import pandas as pd
from rdkit.Chem import AllChem, Descriptors
import boto3
from minio import Minio
from io import BytesIO
import asyncio
from src.molecule.dao import MoleculeDao

# MinIO client setup (to simulate S3)
minio_client = Minio(
    endpoint="storage:9000",  # MinIO service address
    access_key="minio_access_key",
    secret_key="minio_secret_key",
    secure=False  # MinIO is running on HTTP
)

bucket_name = "my_bucket"

# Default arguments for the DAG
default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'start_date': datetime(2023, 10, 17),
    'retries': 1,
}

# Define the DAG
with DAG('molecule_processing_dag', default_args=default_args, schedule_interval='@daily') as dag:

    # Task 1: Extract Data
    def extract_data(**kwargs):
        today = datetime.today().date()
        molecules = asyncio.run(MoleculeDao.get_all_mols())
        molecules_dicts = [
            {"identifier": mol.identifier, "name": mol.name, "smiles": mol.smiles, "created_at": mol.created_at}
            for mol in molecules
        ]
        return molecules_dicts

    # Task 2: Transform Data
    def transform_data(**kwargs):
        ti = kwargs['ti']
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

    # Task 3: Save Data and Upload to MinIO
    def save_data(**kwargs):
        ti = kwargs['ti']
        mols_df = ti.xcom_pull(task_ids='transform_data')

        # Save to my Excel in memory
        output = BytesIO()
        mols_df.to_excel(output, index=False)
        output.seek(0)

        # to upload the file to MinIO
        file_name = f"molecule_properties_{datetime.now().strftime('%Y-%m-%d')}.xlsx"
