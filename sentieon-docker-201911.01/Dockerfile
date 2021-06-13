From ubuntu:18.04

RUN apt-get update && apt-get install -y python \
python-pip \
dnsutils

ENV PATH="/opt/sentieon-genomics-201911.01/bin/:${PATH}"
COPY sentieon-genomics-201911.01 /opt/sentieon-genomics-201911.01
COPY gen_credentials.py /opt/gen_credentials.py
COPY sentieonStartup.sh /opt/sentieonStartup.sh
COPY requirements.txt .
COPY mito-workflow.sh .

RUN pip install -r requirements.txt

CMD ["/bin/bash", "-c", "source mito-workflow.txt"]
