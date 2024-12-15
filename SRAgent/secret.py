# import
## batteries
import os
## 3rd party
from google.cloud import secretmanager

def get_secret(secret_name: str, add_suffix: bool=True) -> str:
    """
    Fetch secret from GCP Secret Manager.
    Rquired environment variables: GCP_PROJECT_ID, GOOGLE_APPLICATION_CREDENTIALS
    Args:
        secret_id: The secret id
        add_suffix: Add suffix to secret name based on PY_CONFIG_ACTIVE
    Returns:
        The secret value
    """
    # get secret version
    if add_suffix:
        secret_name = add_secret_suffix(secret_name)
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
    suffix = "TEST" if os.getenv("PY_CONFIG_ACTIVE") == "TEST" else "PROD"
    return f"{secret_name}_{suffix}"

def get_env(env_var: str, add_suffix: bool=True) -> str:
    """
    """
    if os.getenv(env_var):
        return os.environ[env_var]
    if add_suffix:
        suffix = "TEST" if os.getenv("PY_CONFIG_ACTIVE") == "TEST" else "PROD"
        env_var = f"{env_var}_{suffix}"
    return os.environ[env_var]

# main
if __name__ == '__main__':
    print(get_secret(f"GCP_SQL_DB_PASSWORD"))