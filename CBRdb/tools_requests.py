from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util import Retry


def prepare_session():
    """
    Prepares a requests session with retry logic.

    Returns:
    requests.Session: A configured session object with retry logic.
    """
    # Make the session
    s = Session()
    # Add retries
    retries = Retry(
        total=5,
        backoff_factor=0.1,
        status_forcelist=[502, 503, 504],
        allowed_methods={'POST'},
    )
    # Mount the session
    s.mount('https://', HTTPAdapter(max_retries=retries))
    return s
