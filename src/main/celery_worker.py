from celery import Celery

celery = Celery(
    'main.celery_worker',
    broker='redis://redis:6379/0',
    backend='redis://redis:6379/0'
)

celery.conf.update(task_track_started=True)
celery.autodiscover_tasks(['main.tasks'])
