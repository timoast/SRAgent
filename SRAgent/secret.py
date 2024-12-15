# import
## batteries
import os
## 3rd party
from google.cloud import secretmanager

def get_secret(secret_name: str, add_suffix: bool=True, return_env: bool=True) -> str:
    """
    Fetch secret from GCP Secret Manager.
    Rquired environment variables: GCP_PROJECT_ID, GOOGLE_APPLICATION_CREDENTIALS
    Args:
        secret_id: The secret id
        add_suffix: Add suffix to secret name based on PY_CONFIG_ACTIVE
        return_env: Return secret value from env if already present
    Returns:
        The secret value
    """
    # if secret is already in env, return it
    if return_env and os.getenv(secret_name):
        return os.environ[secret_name]
    # get secret version
    if add_suffix:
        secret_name = add_secret_suffix(secret_name)
        if return_env and os.getenv(secret_name):
            return os.environ[secret_name]
    # fetch secret from GCP Secret Manager
    client = secretmanager.SecretManagerServiceClient()
    name = f"projects/{os.environ['GCP_PROJECT_ID']}/secrets/{secret_name}/versions/latest"
    response = client.access_secret_version(request={"name": name})
    return response.payload.data.decode('UTF-8')

def add_secret_suffix(secret_name: str) -> str:
    """
    Add suffix to secret name based on active config.
    Args:
        secret_name: The secret name
    Returns:
        The secret name with suffix
    """
    suffix = "TEST" if os.getenv("PY_CONFIG_ACTIVE", "").upper() == "TEST" else "PROD"
    return f"{secret_name}_{suffix}"

# main
if __name__ == '__main__':
    print(get_secret(f"GCP_SQL_DB_PASSWORD"))