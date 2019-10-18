from django.http import HttpResponse
from django.template import loader
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.core.mail import send_mail, BadHeaderError
from django.conf import settings

#@csrf_exempt
def home(request):
	template = loader.get_template('home/home.html')
	return HttpResponse(template.render())

@csrf_exempt
def mail(request):
	try:
		name = request.POST['Name']
		email = request.POST['Email']
		subject = request.POST['Subject']
		message = request.POST['Message']
	except KeyError:
		return HttpResponse('Mail could not be sent')

	send_mail(subject, message, settings.EMAIL_HOST_USER, ['michiel.camps@gmail.com'])

	return render(request, 'home/mailus.html', {'name': name, 'email': email,'subject': subject, 'message': message})
