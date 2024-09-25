from celery import Celery
import logging

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    filename="/app/celery_logs.log",  # Use an absolute path within the Docker container
    filemode="w",
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# Create a logger and ensure it doesn't propagate to the root logger
logger = logging.getLogger(__name__)
logger.propagate = False

# Initialize Celery
celery = Celery(
    'tasks',
    broker='redis://redis:6379/0',
    backend='redis://redis:6379/0'
)

celery.conf.update(task_track_started=True)
celery.autodiscover_tasks(['src.molecule.celery_tasks'])
