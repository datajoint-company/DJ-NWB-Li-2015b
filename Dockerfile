FROM datajoint/jupyter:python3.6

RUN pip uninstall -y datajoint
RUN pip install "git+https://github.com/datajoint/datajoint-python.git@dev#egg=datajoint"
RUN pip install pynwb
ADD . /src/li2015b
RUN pip install -e /src/li2015b
