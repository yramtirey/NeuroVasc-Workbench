import { defineConfig } from '@playwright/test';
export default defineConfig({
  testDir: './e2e', timeout: 60000, workers: 1,
  use: { baseURL: 'http://localhost:5173', headless: true, launchOptions: { args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'] } },
  webServer: [
    { command: '../.venv/bin/python -m uvicorn backend.app.main:app --app-dir .. --host 127.0.0.1 --port 8000', url: 'http://127.0.0.1:8000', timeout: 180000, reuseExistingServer: !process.env.CI },
    { command: 'npm run dev -- --host localhost', url: 'http://localhost:5173', timeout: 180000, reuseExistingServer: !process.env.CI },
  ],
});
