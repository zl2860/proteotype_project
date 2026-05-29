#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse


HTML = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>ST10 Training Monitor</title>
  <style>
    :root {
      color-scheme: light;
      --bg: #f7f7f4;
      --ink: #202322;
      --muted: #666b68;
      --line: #d8ddd7;
      --panel: #ffffff;
      --accent: #006c67;
      --accent2: #8c3f2a;
      --accent3: #576f9e;
      --warn: #b45b00;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: var(--bg);
      color: var(--ink);
    }
    header {
      display: flex;
      align-items: center;
      gap: 16px;
      padding: 16px 20px;
      border-bottom: 1px solid var(--line);
      background: var(--panel);
      position: sticky;
      top: 0;
      z-index: 2;
    }
    h1 {
      margin: 0;
      font-size: 18px;
      font-weight: 650;
      white-space: nowrap;
    }
    select, button {
      height: 34px;
      border: 1px solid var(--line);
      background: #fff;
      color: var(--ink);
      border-radius: 6px;
      padding: 0 10px;
      font: inherit;
    }
    select { min-width: min(720px, 55vw); }
    main {
      padding: 18px 20px 28px;
      display: grid;
      gap: 18px;
    }
    .metrics {
      display: grid;
      grid-template-columns: repeat(5, minmax(120px, 1fr));
      gap: 10px;
    }
    .metric {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 12px;
      min-height: 76px;
    }
    .metric .label {
      color: var(--muted);
      font-size: 12px;
      margin-bottom: 8px;
    }
    .metric .value {
      font-size: 22px;
      font-weight: 700;
      font-variant-numeric: tabular-nums;
    }
    .grid {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 18px;
    }
    section {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
      min-width: 0;
    }
    h2 {
      font-size: 14px;
      margin: 0 0 10px;
      color: var(--muted);
      font-weight: 650;
    }
    canvas {
      width: 100%;
      height: 280px;
      display: block;
    }
    table {
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
      background: var(--panel);
    }
    th, td {
      border-bottom: 1px solid var(--line);
      padding: 7px 8px;
      text-align: right;
      font-variant-numeric: tabular-nums;
    }
    th:first-child, td:first-child { text-align: left; }
    .status {
      color: var(--muted);
      margin-left: auto;
      font-size: 13px;
      white-space: nowrap;
    }
    .legend {
      display: flex;
      gap: 14px;
      flex-wrap: wrap;
      color: var(--muted);
      font-size: 12px;
      margin-top: 8px;
    }
    .dot {
      display: inline-block;
      width: 10px;
      height: 10px;
      border-radius: 50%;
      margin-right: 5px;
      vertical-align: -1px;
    }
    @media (max-width: 980px) {
      header { align-items: stretch; flex-direction: column; }
      select { min-width: 0; width: 100%; }
      .status { margin-left: 0; }
      .metrics { grid-template-columns: repeat(2, minmax(120px, 1fr)); }
      .grid { grid-template-columns: 1fr; }
    }
  </style>
</head>
<body>
  <header>
    <h1>ST10 Training Monitor</h1>
    <select id="runSelect" aria-label="Run"></select>
    <button id="refreshButton">Refresh</button>
    <div class="status" id="status">Loading...</div>
  </header>
  <main>
    <div class="metrics">
      <div class="metric"><div class="label">Epochs</div><div class="value" id="epochs">-</div></div>
      <div class="metric"><div class="label">Best Val Loss</div><div class="value" id="bestVal">-</div></div>
      <div class="metric"><div class="label">Latest Val LM</div><div class="value" id="valLm">-</div></div>
      <div class="metric"><div class="label">Latest Val Time</div><div class="value" id="valTime">-</div></div>
      <div class="metric"><div class="label">Train-Val Gap</div><div class="value" id="gap">-</div></div>
    </div>
    <div class="grid">
      <section>
        <h2>Total Loss</h2>
        <canvas id="lossCanvas" width="900" height="320"></canvas>
        <div class="legend">
          <span><i class="dot" style="background:var(--accent)"></i>train_loss</span>
          <span><i class="dot" style="background:var(--accent2)"></i>val_loss</span>
        </div>
      </section>
      <section>
        <h2>Trajectory Losses</h2>
        <canvas id="taskCanvas" width="900" height="320"></canvas>
        <div class="legend">
          <span><i class="dot" style="background:var(--accent3)"></i>val_lm_loss</span>
          <span><i class="dot" style="background:var(--warn)"></i>val_time_loss</span>
        </div>
      </section>
    </div>
    <section>
      <h2>History</h2>
      <table>
        <thead id="thead"></thead>
        <tbody id="tbody"></tbody>
      </table>
    </section>
  </main>
  <script>
    const runSelect = document.getElementById('runSelect');
    const statusEl = document.getElementById('status');
    const fields = ['epoch','train_loss','val_loss','train_lm_loss','val_lm_loss','train_time_loss','val_time_loss'];
    let timer = null;

    function fmt(x) {
      if (x === undefined || x === null || Number.isNaN(Number(x))) return '-';
      return Number(x).toFixed(4);
    }

    async function fetchJson(url) {
      const r = await fetch(url, {cache: 'no-store'});
      if (!r.ok) throw new Error(await r.text());
      return await r.json();
    }

    async function loadRuns() {
      const previous = runSelect.value;
      const runs = await fetchJson('/api/runs');
      runSelect.innerHTML = '';
      for (const run of runs) {
        const opt = document.createElement('option');
        opt.value = run.history;
        opt.textContent = `${run.name} (${run.epochs} epochs)`;
        runSelect.appendChild(opt);
      }
      if (runs.length) runSelect.value = runs.some(r => r.history === previous) ? previous : runs[0].history;
      await loadHistory();
    }

    async function loadHistory() {
      if (!runSelect.value) return;
      const data = await fetchJson(`/api/history?path=${encodeURIComponent(runSelect.value)}`);
      render(data.rows || []);
      statusEl.textContent = `Updated ${new Date().toLocaleTimeString()} | ${data.path}`;
    }

    function setMetric(id, value) { document.getElementById(id).textContent = value; }

    function render(rows) {
      const latest = rows[rows.length - 1] || {};
      const best = rows.reduce((m, r) => Math.min(m, Number(r.val_loss ?? Infinity)), Infinity);
      setMetric('epochs', rows.length || '-');
      setMetric('bestVal', Number.isFinite(best) ? fmt(best) : '-');
      setMetric('valLm', fmt(latest.val_lm_loss));
      setMetric('valTime', fmt(latest.val_time_loss));
      const gap = latest.val_loss !== undefined && latest.train_loss !== undefined ? latest.val_loss - latest.train_loss : null;
      setMetric('gap', gap === null ? '-' : fmt(gap));
      drawChart(document.getElementById('lossCanvas'), rows, [
        ['train_loss', getCss('--accent')],
        ['val_loss', getCss('--accent2')],
      ]);
      drawChart(document.getElementById('taskCanvas'), rows, [
        ['val_lm_loss', getCss('--accent3')],
        ['val_time_loss', getCss('--warn')],
      ]);
      renderTable(rows);
    }

    function getCss(name) {
      return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    }

    function drawChart(canvas, rows, series) {
      const ctx = canvas.getContext('2d');
      const w = canvas.width, h = canvas.height;
      ctx.clearRect(0, 0, w, h);
      ctx.fillStyle = '#fff';
      ctx.fillRect(0, 0, w, h);
      const pad = {l: 54, r: 18, t: 18, b: 38};
      const values = [];
      for (const [key] of series) for (const r of rows) if (r[key] !== undefined) values.push(Number(r[key]));
      if (!rows.length || !values.length) return;
      let minY = Math.min(...values), maxY = Math.max(...values);
      const span = Math.max(maxY - minY, 1e-6);
      minY -= span * 0.08; maxY += span * 0.08;
      const xFor = i => pad.l + (rows.length === 1 ? 0 : i * (w - pad.l - pad.r) / (rows.length - 1));
      const yFor = v => h - pad.b - (Number(v) - minY) * (h - pad.t - pad.b) / (maxY - minY);
      ctx.strokeStyle = '#d8ddd7';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, h - pad.b); ctx.lineTo(w - pad.r, h - pad.b);
      ctx.stroke();
      ctx.fillStyle = '#666b68';
      ctx.font = '12px system-ui';
      for (let i = 0; i <= 4; i++) {
        const y = pad.t + i * (h - pad.t - pad.b) / 4;
        const v = maxY - i * (maxY - minY) / 4;
        ctx.fillText(v.toFixed(2), 8, y + 4);
        ctx.strokeStyle = '#eef0ed';
        ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();
      }
      for (const [key, color] of series) {
        ctx.strokeStyle = color; ctx.lineWidth = 2.5; ctx.beginPath();
        rows.forEach((r, i) => {
          const x = xFor(i), y = yFor(r[key]);
          if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        });
        ctx.stroke();
        ctx.fillStyle = color;
        rows.forEach((r, i) => {
          const x = xFor(i), y = yFor(r[key]);
          ctx.beginPath(); ctx.arc(x, y, 3, 0, Math.PI * 2); ctx.fill();
        });
      }
      ctx.fillStyle = '#666b68';
      rows.forEach((r, i) => {
        if (i === 0 || i === rows.length - 1 || rows.length <= 10) ctx.fillText(r.epoch, xFor(i) - 4, h - 12);
      });
    }

    function renderTable(rows) {
      document.getElementById('thead').innerHTML = `<tr>${fields.map(f => `<th>${f}</th>`).join('')}</tr>`;
      document.getElementById('tbody').innerHTML = rows.slice().reverse().map(r => (
        `<tr>${fields.map(f => `<td>${f === 'epoch' ? (r[f] ?? '-') : fmt(r[f])}</td>`).join('')}</tr>`
      )).join('');
    }

    document.getElementById('refreshButton').addEventListener('click', loadHistory);
    runSelect.addEventListener('change', loadHistory);
    loadRuns().catch(err => statusEl.textContent = err.message);
    timer = setInterval(() => loadHistory().catch(err => statusEl.textContent = err.message), 3000);
    setInterval(() => loadRuns().catch(err => statusEl.textContent = err.message), 30000);
  </script>
</body>
</html>
"""


class DashboardHandler(BaseHTTPRequestHandler):
    root: Path = Path(".")

    def json_response(self, data, status: int = 200) -> None:
        payload = json.dumps(data, ensure_ascii=False).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def text_response(self, text: str, status: int = 200, content_type: str = "text/html; charset=utf-8") -> None:
        payload = text.encode()
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/":
            self.text_response(HTML)
        elif parsed.path == "/api/runs":
            self.json_response(self.list_runs())
        elif parsed.path == "/api/history":
            query = parse_qs(parsed.query)
            path = query.get("path", [""])[0]
            self.json_response(self.read_history(path))
        else:
            self.text_response("Not found", status=404, content_type="text/plain; charset=utf-8")

    def list_runs(self):
        runs = []
        for hist in sorted((self.root / "model_runs").glob("*/history.jsonl"), key=lambda p: p.stat().st_mtime, reverse=True):
            rows = self.parse_history(hist)
            runs.append(
                {
                    "name": hist.parent.name,
                    "history": str(hist.relative_to(self.root)),
                    "epochs": len(rows),
                    "mtime": hist.stat().st_mtime,
                }
            )
        return runs

    def read_history(self, rel_path: str):
        path = (self.root / rel_path).resolve()
        root = self.root.resolve()
        if root not in path.parents:
            return {"path": rel_path, "rows": [], "error": "path outside workspace"}
        if not path.exists():
            return {"path": rel_path, "rows": [], "error": "history file not found"}
        return {"path": rel_path, "rows": self.parse_history(path)}

    @staticmethod
    def parse_history(path: Path):
        rows = []
        with path.open() as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    rows.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
        return rows

    def log_message(self, fmt, *args):
        return


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    DashboardHandler.root = args.root.resolve()
    server = ThreadingHTTPServer((args.host, args.port), DashboardHandler)
    print(f"Training dashboard: http://{args.host}:{args.port}", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
