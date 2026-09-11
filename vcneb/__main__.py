"""Allow `python -m vcneb` to use the package CLI."""

from .cli import main


if __name__ == "__main__":
    raise SystemExit(main())
