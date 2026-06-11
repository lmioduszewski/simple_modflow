"""Shared display helpers used by the canonical master notebook set."""

from IPython.display import HTML, display


CANONICAL_CSS = """
<style>
:root { --ink:#18332f; --pine:#276257; --water:#2f7fa3; --sand:#f2e5c4; --mist:#eef5f2; }
.canonical-hero {
  padding: 28px 32px; border-radius: 18px; color: white;
  background: linear-gradient(125deg, #173d37 0%, #27726a 54%, #3d8ba3 100%);
  box-shadow: 0 12px 30px rgba(23,61,55,.18); margin: 8px 0 24px;
}
.canonical-hero h1 { margin:0 0 8px; font-size:2.1rem; letter-spacing:.01em; }
.canonical-hero p { margin:0; max-width:850px; font-size:1.04rem; line-height:1.55; }
.canonical-card {
  padding: 16px 20px; border-radius: 12px; background:var(--mist);
  border-left:5px solid var(--water); margin:14px 0;
}
.canonical-card strong { color:var(--pine); }
.canonical-grid { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:12px; }
.canonical-chip { padding:9px 12px; border-radius:999px; background:var(--sand); color:var(--ink); }
</style>
"""


def notebook_header(number: str, title: str, subtitle: str) -> None:
    """Render the shared canonical notebook header."""

    display(
        HTML(
            CANONICAL_CSS
            + f"""
            <div class="canonical-hero">
              <h1>{number} · {title}</h1>
              <p>{subtitle}</p>
            </div>
            """
        )
    )
