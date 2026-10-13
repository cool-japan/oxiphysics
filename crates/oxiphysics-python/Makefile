.PHONY: dev test rust-test build clean

dev:
	maturin develop --release

test: dev
	python -m pytest python/tests/ -v

rust-test:
	cargo nextest run -p oxiphysics-python

build:
	maturin build --release --out dist/

clean:
	rm -rf dist/ target/wheels/
