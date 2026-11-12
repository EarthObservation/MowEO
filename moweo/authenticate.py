from .auth import connect_project_account


if __name__ == "__main__":
    connect_project_account(force_interactive=True, confirm=True)
    print("Project-local authentication saved successfully.")
