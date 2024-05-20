#!/bin/bash

#SBATCH --job-name={{ job_name }}
#SBATCH --nodes={{ nodes }}
#SBATCH --ntasks-per-node={{ ntasks_per_node }}
#SBATCH --output={{ output_file }}
#SBATCH --error={{ error_file }}
#SBATCH --time={{ time }}
#SBATCH --cluster={{ cluster }}
#SBATCH --partition={{ partition }}
{% if account %}
#SBATCH --account={{ account }}
{% endif %}

module purge
{% for module in modules %}
module load {{ module }}
{% endfor %}

{% for key, value in env_vars.items() %}
export {{ key }}={{ value }}
{% endfor %}

{% for command in commands_pre %}
command
{% endfor %}

{% for command in commands_run %}
command
{% endfor %}

{% for command in commands_post %}
command
{% endfor %}
