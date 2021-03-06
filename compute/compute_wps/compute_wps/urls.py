from __future__ import absolute_import

from django.conf import settings
from django.conf.urls import include
from django.urls import re_path
from rest_framework import routers
from rest_framework.schemas import get_schema_view
from rest_framework.urlpatterns import format_suffix_patterns

from compute_wps.views import auth
from compute_wps.views import job
from compute_wps.views import metrics
from compute_wps.views import service

router = routers.SimpleRouter()
router.register("message", job.MessageViewSet)
router.register("process", job.ProcessViewSet)

status_router = routers.SimpleRouter()
status_router.register("status", job.StatusViewSet)

output_router = routers.SimpleRouter()
output_router.register("output", job.OutputViewSet)

job_router = routers.SimpleRouter()
job_router.register("job", job.JobViewSet)

# Only require the format on job-detail view
router_urls = router.urls
router_urls.append(job_router.urls[0])
router_urls.append(status_router.urls[0])
router_urls.append(output_router.urls[0])
router_urls.extend(
    format_suffix_patterns(job_router.urls[1:2], allowed=["json", "wps"])
)
router_urls.extend(
    format_suffix_patterns(status_router.urls[1:2], allowed=["json", "wps"])
)

schema = get_schema_view(
    title="WPS API", patterns=router.urls, url=f"/wps/api"
    # title="WPS API", patterns=router.urls, url=f"{settings.BASE_URL}/wps/api"
)


api_urlpatterns = [
    re_path("^", include(router_urls)),
    re_path("^schema/?$", schema),
    re_path("^ping/?$", service.ping),
    re_path("^metrics/?$", metrics.metrics_view),
]

auth_urlpatterns = [
    re_path("^client_registration/?$", auth.client_registration),
    re_path("^login/?$", auth.login),
    re_path("^oauth_callback/?$", auth.login_complete),
]

wps_urlpatterns = [
    re_path("^$", service.wps_entrypoint),
    re_path("^api/", include(api_urlpatterns)),
    re_path("^auth/", include(auth_urlpatterns)),
]

urlpatterns = [
    re_path("^wps/?", include(wps_urlpatterns)),
]

if settings.DEBUG:
    from django.conf.urls.static import static

    urlpatterns += static(
        settings.STATIC_URL, document_root=settings.STATIC_ROOT
    )
