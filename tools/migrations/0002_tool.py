# Generated by Django 2.2.6 on 2019-10-11 06:54

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('tools', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Tool',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('tool_name', models.CharField(max_length=30)),
                ('tool_date', models.DateTimeField(verbose_name='Published online')),
                ('tool_description', models.CharField(max_length=1000)),
                ('tool_user_number', models.IntegerField(default=0)),
            ],
        ),
    ]