import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

import { nodePolyfills } from "@blocksquaredev/vite-plugin-node-polyfills";


export default defineConfig({
  server: {
    port: 5173,
    strictPort: true,
  },
  plugins: [
    react(),

    nodePolyfills({
      include: [
        "events",
      ],

      globals: {
        Buffer: true,
        global: true,
        process: true,
      },

      protocolImports: true,
    }),
  ],
});
