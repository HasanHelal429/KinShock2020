---
description: Summarize project state — runs, recent findings, anything running
---
Give me the current project state:
1. `git log --oneline -8` and `git status --short`.
2. List `runs/` and which have completed diags.
3. The last ~40 lines of RESULTS.md (most recent findings).
4. Any WarpX process currently running (`pgrep -af warpx.1d`) and, if a progress.log
   exists in its run dir, its latest checkpoint line.
Summarize concisely: what's done, what's in flight, what's next.
