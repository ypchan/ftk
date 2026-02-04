import typer
from ftk.commands import stat

app = typer.Typer(
    name="ftk",
    help="ftk: extensible CLI toolkit.",
    add_completion=False,
)

# Register subcommands here
app.add_typer(stat.app, name="stat")

def main():
    app()

if __name__ == "__main__":
    main()
