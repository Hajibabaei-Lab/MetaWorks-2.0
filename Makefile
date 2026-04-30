.PHONY: dev dev-backend dev-frontend install install-frontend test test-backend test-frontend lint lint-backend lint-frontend build

dev-backend:
	uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload

dev-frontend:
	cd frontend && npm run dev

dev:
	@echo "Starting backend and frontend..."
	@$(MAKE) dev-backend & $(MAKE) dev-frontend

install-frontend:
	cd frontend && npm ci

install: install-frontend

test-backend:
	PYTHONPATH=$(PWD) pytest -v

test-frontend:
	cd frontend && npm run test

test: test-backend test-frontend

lint-backend:
	ruff check .

lint-frontend:
	cd frontend && npx vue-tsc --noEmit 2>/dev/null || true

lint: lint-backend lint-frontend

build:
	cd frontend && npm run build
