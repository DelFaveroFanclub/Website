from django.urls import path

from . import views

urlpatterns = [
	path('',views.index,name='index'),
    path('<str:tool_name>/', views.detail, name='detail'),
    path('<str:tool_name>/results/', views.results,name='results'),
    path('<str:tool_name>/upload/', views.upload, name='upload'),
]
