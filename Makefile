init:
	pip install -r requirements.txt
	python setup.py build_ext --inplace

test:
	nosetests tests