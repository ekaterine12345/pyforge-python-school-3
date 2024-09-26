import boto3
import json

# Initialize a session
session = boto3.Session()

# Initialize the Lambda client from session
lambda_client = session.client('lambda', region_name='us-east-2')

# define the event
event = {
    "names": ["Kate", "Katerine", "Alice", "Bob", "James"]
}

# Invoke the Lambda function
try:
    response = lambda_client.invoke(
        FunctionName='HelloStudentFunction',
        InvocationType='RequestResponse',  # RequestResponse waits for the function to complete
        Payload=json.dumps(event),
    )

    # Read the response payload
    response_payload = json.loads(response['Payload'].read())
    print("Lambda Response:", response_payload)

except Exception as e:
    print(f"Error invoking Lambda function: {e}")
