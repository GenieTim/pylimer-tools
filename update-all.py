import warnings
import random
from subprocess import call

import pkg_resources

packages = [dist.project_name for dist in pkg_resources.working_set]
random.shuffle(packages)
call("pip install --upgrade " + ' '.join(packages), shell=True)
# for package in packages:
#     try:
#         call("pip install --upgrade " + package, shell=True)
#     except Exception as e:
#         warnings.warn(e)
