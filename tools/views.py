from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.urls import reverse
from .models import Tool
import os
from locomotif.settings import MEDIA_ROOT
from .LocoMotif import LocoMotif
from .LocoMotif import Gene

def index(request):
    tool_list = Tool.objects.all()
    context = {'tool_list': tool_list}
    return render(request, 'tools/index.html', context)

def detail(request,tool_name):
    tool = get_object_or_404(Tool,pk=Tool.objects.get(tool_name=tool_name).id)
    return render(request,'tools/detail.html',{'tool':tool})


def results(request, tool_name):
    tool = get_object_or_404(Tool, pk=Tool.objects.get(tool_name=tool_name).id)
    filename = tool.file.name
    full_path = os.path.join(MEDIA_ROOT, filename)
    freq, expfreq = LocoMotif.FindMotif(file=full_path,motif=Gene.Motif(seq=tool.motif),overlap=float(tool.overlap))
    return render(request, 'tools/results.html',{'tool': tool, 'freq': freq*100, 'expfreq': expfreq*100})

def upload(request, tool_name):
    tool = get_object_or_404(Tool,pk=Tool.objects.get(tool_name=tool_name).id)
    tool.file = request.FILES["file"]
    tool.motif = request.POST['motif']
    tool.overlap = request.POST['overlap']
    tool.save()
    return HttpResponseRedirect(reverse('results', args=(tool.tool_name,)))
