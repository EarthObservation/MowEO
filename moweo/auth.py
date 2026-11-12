from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

from .config import AUTH_ACCOUNT_FILE, AUTH_DIR, AUTH_TOKEN_FILE, OPENEO_BACKEND_URL, YEAR


class ProjectRefreshTokenStore:
    """Minimal refresh-token store compatible with openEO Connection.

    Tokens are intentionally isolated inside the configured year workspace instead
    of the global user profile. This allows concurrent runs with different CDSE
    accounts in separate projects/years.
    """

    def __init__(self, path: Path):
        self.path = Path(path)

    def _read(self) -> dict:
        if not self.path.exists():
            return {}
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
            return data if isinstance(data, dict) else {}
        except Exception:
            return {}

    def _write(self, data: dict) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(self.path.suffix + ".part")
        tmp.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
        tmp.replace(self.path)
        try:
            os.chmod(self.path, 0o600)
        except Exception:
            pass

    def get_refresh_token(self, issuer: str, client_id: str):
        item = self._read().get(str(issuer), {}).get(str(client_id), {})
        if isinstance(item, dict):
            return item.get("refresh_token")
        return None

    def set_refresh_token(self, issuer: str, client_id: str, refresh_token: str) -> None:
        data = self._read()
        data.setdefault(str(issuer), {})[str(client_id)] = {
            "date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "refresh_token": str(refresh_token),
        }
        self._write(data)

    def remove(self) -> None:
        if self.path.exists():
            self.path.unlink()


def connect_project_account(*, force_interactive: bool = True, confirm: bool = True):
    """Authenticate against openEO with a project-local token store.

    By default the device login is forced at the beginning of every acquisition
    invocation. The resulting refresh token is stored only in the current year
    workspace and can be used by the same long-running Connection for automatic
    token renewal. It never touches the default global openEO token cache.
    """
    from openeo.rest.connection import Connection

    AUTH_DIR.mkdir(parents=True, exist_ok=True)
    store = ProjectRefreshTokenStore(AUTH_TOKEN_FILE)
    connection = Connection(OPENEO_BACKEND_URL, refresh_token_store=store)

    if force_interactive:
        print("Interactive CDSE/openEO login is required for this acquisition run.")
        print("This project uses its own token store:", AUTH_TOKEN_FILE)
        connection = connection.authenticate_oidc_device(store_refresh_token=True)
    else:
        connection = connection.authenticate_oidc(store_refresh_token=True)

    account = connection.describe_account()
    user_id = str(account.get("user_id", "")).strip()
    name = str(account.get("name", "") or "").strip()
    if not user_id:
        raise RuntimeError(f"Could not read openEO user_id from account: {account}")

    print("Authentication: OK")
    print("Year         :", YEAR)
    print("openEO user  :", user_id)
    if name:
        print("Account name :", name)

    if confirm:
        answer = input("Type YES to use this account for this acquisition run: ").strip()
        if answer != "YES":
            raise SystemExit("Account not confirmed. No acquisition jobs were submitted.")

    AUTH_ACCOUNT_FILE.write_text(
        json.dumps(
            {
                "year": YEAR,
                "user_id": user_id,
                "name": name,
                "backend": OPENEO_BACKEND_URL,
                "confirmed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return connection
