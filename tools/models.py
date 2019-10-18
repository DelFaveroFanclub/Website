import datetime

from django.db import models
from django.utils import timezone
#from django import forms

class Question(models.Model):
	question_text = models.CharField(max_length=200)
	pub_date = models.DateTimeField('date published')

	def __str__(self):
		return self.question_text

	def was_published_recently(self):
		return self.pub_date >= timezone.now() - datetime.timedelta(days=1)

class Choice(models.Model):
	question = models.ForeignKey(Question,on_delete=models.CASCADE)
	choice_text = models.CharField(max_length=200)
	votes=models.IntegerField(default=0)

	def __str__(self):
		return self.choice_text

class Tool(models.Model):
	tool_name = models.CharField(max_length=30)
	tool_date = models.DateTimeField('Published online')
	tool_description = models.CharField(max_length=1000)
	tool_user_number = models.IntegerField(default=0)
	file = models.FileField(default="",upload_to='uploads/')
	motif = models.CharField(max_length=30,default='')
	overlap = models.CharField(max_length=10,default='')

	def __str__(self):
		return self.tool_name
