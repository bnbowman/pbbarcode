.PHONY: doc doc-clean

SHELL = /bin/bash -e

all: build install

build:
	python setup.py build --executable="/usr/bin/env python"

bdist:
	python setup.py build --executable="/usr/bin/env python"
	python setup.py bdist --formats=egg

install:
	python setup.py install

develop:
	python setup.py develop

test:
	find tests -name "*.py" | xargs nosetests
	find tests/cram -name "*.t" | grep -v consensus.t | xargs cram --verbose 

clean: doc-clean
	rm -rf build/;\
	find . -name "*.egg-info" | xargs rm -rf;\
	find . -name "*.pyc" | xargs rm -rf;\
	rm -rf dist/
	make -C src/C clean

doc-clean:
	make -C doc clean

doc:
	make -C doc html

pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'pbtools.barcode=='>/dev/null \
      && pip uninstall -y pbtools.barcode \
      || echo -n ''
	@pip freeze|grep 'pbbarcode=='>/dev/null \
      && pip uninstall -y pbbarcode \
      || echo -n ''
	@pip install --no-index \
          --install-option="--install-scripts=$(PREFIX)/bin" \
          ./
