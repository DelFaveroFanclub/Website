# Generated by Django 2.2.6 on 2019-10-14 13:29

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('tools', '0006_tool_file'),
    ]

    operations = [
        migrations.AddField(
            model_name='tool',
            name='motif',
            field=models.CharField(default='', max_length=30),
        ),
        migrations.AddField(
            model_name='tool',
            name='overlap',
            field=models.CharField(default='', max_length=10),
        ),
    ]