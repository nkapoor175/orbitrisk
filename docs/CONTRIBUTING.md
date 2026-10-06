# Git workflow

Navika is the maintainer. Rules, kept short on purpose:

1. `main` always runs. Nobody commits to it directly.
2. Work on your own branch: `git checkout -b parthvi/monte-carlo`, `vishi/error-model`, `navika/download`.
3. Commit small and often with a clear message: `Add residual builder`, not `update`.
4. When a script works, push your branch and open a pull request. Navika reviews and merges.
5. Never commit `.env`, passwords, or anything in `data/raw/`.
6. Only edit scripts you own (see the table in the README). Need a change in someone else's script? Ask them or open an issue.
7. Hand-off CSVs in `data/processed/` must follow `docs/DATA_CONTRACT.md` exactly.

First-time clone:

```bash
git clone <repo-url>
cd orbitrisk
git checkout -b yourname/first-task
```
