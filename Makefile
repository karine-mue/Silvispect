.PHONY: help test lint format check demo data clean

help:
	@echo "make test    - run the test suite"
	@echo "make lint    - ruff check and format check"
	@echo "make format  - apply ruff formatting"
	@echo "make check   - lint and test"
	@echo "make demo    - generate a stand and inspect it"
	@echo "make data    - regenerate the sample plot in data/"
	@echo "make clean   - remove caches and build artefacts"

test:
	pytest

lint:
	ruff check .
	ruff format --check .
	mypy silvispect

format:
	ruff format .
	ruff check --fix .

check: lint test

demo:
	python -m silvispect demo

data:
	python scripts/make_sample_data.py

clean:
	rm -rf .pytest_cache .ruff_cache build dist *.egg-info
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
