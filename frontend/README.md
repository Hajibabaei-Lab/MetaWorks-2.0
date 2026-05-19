# MetaWorks UI

Standalone Vue 3 frontend for the MetaWorks FastAPI backend.

## Stack

- Vue 3 + Vite + TypeScript
- Vue Router
- TanStack Query
- Zod
- Element Plus
- Vitest + Playwright
- Caddy for production SPA serving and `/api` reverse proxying

## Development

1. Start the backend from `MetaWorks-2.0` on `http://127.0.0.1:8000`
2. Install dependencies:
   ```bash
   npm install
   ```
3. Start the frontend:
   ```bash
   npm run dev
   ```
4. Open `http://127.0.0.1:5173`

The Vite dev server proxies `/api/*` to the backend.

## Type Generation

Regenerate API types from the backend OpenAPI schema:

```bash
npm run generate:types
```

## Production Container

Build and run the UI container:

```bash
docker build -t metaworks-ui .
docker run --rm -p 8080:8080 metaworks-ui
```

The container expects a Compose network alias named `backend` for API proxying.

