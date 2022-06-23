FROM public.ecr.aws/lambda/python:3.9

COPY ./ ${LAMBDA_TASK_ROOT}

COPY requirements.txt .


RUN yum install -y yum-utils rpmdevtools

WORKDIR /tmp

RUN yumdownloader --resolve expat glib2 libffi libffi-devel cairo pango && rpmdev-extract *rpm

WORKDIR ${LAMBDA_TASK_ROOT}

RUN cp -P -R /tmp/*/usr/lib64/* ${LAMBDA_TASK_ROOT}

RUN ln libgobject-2.0.so.0 libgobject-2.0.so && \
    ln libcairo.so.2 libcairo.so && \
    ln libpango-1.0.so.0 pango-1.0 && \
    ln libpangoft2-1.0.so.0 pangoft2-1.0 && \
    ln libpangocairo-1.0.so.0 pangocairo-1.0

RUN pip3 install -r requirements.txt --target "${LAMBDA_TASK_ROOT}"


CMD [ "app.lambda_handler" ]