from setuptools import setup
import os


def find_stubs(package):
    stubs = []
    for root, dirs, files in os.walk(package):
        for file in files:
            path = os.path.join(root, file).replace(package + os.sep, '', 1)
            stubs.append(path)
    return dict(package=stubs)


setup(
    name='pylimer_tools_cpp-stubs',
    maintainer="pylimer_tools_cpp Developers",
    maintainer_email="example@python.org",
    description="PEP 561 type stubs for pylimer_tools_cpp",
    version='1.0',
    packages=['pylimer_tools_cpp-stubs'],
    # PEP 561 requires these
    install_requires=['pylimer_tools_cpp'],
    package_data=find_stubs('pylimer_tools_cpp-stubs'),
)