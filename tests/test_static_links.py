from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlparse


ROOT = Path(__file__).resolve().parents[1]


class AssetParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.references = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "a" and attrs.get("href"):
            self.references.append(("href", attrs["href"]))
        if tag == "img" and attrs.get("src"):
            self.references.append(("src", attrs["src"]))


def local_references(html):
    parser = AssetParser()
    parser.feed(html)
    for attr, value in parser.references:
        parsed = urlparse(value)
        if parsed.scheme or value.startswith(("#", "mailto:", "tel:")):
            continue
        yield attr, parsed.path


def test_landing_page_references_existing_local_assets():
    html = (ROOT / "index.html").read_text(encoding="utf-8")

    assert "{{" not in html
    assert "C:\\" not in html
    assert "D:\\" not in html

    missing = []
    for attr, ref in local_references(html):
        target = (ROOT / ref).resolve()
        if not target.exists():
            missing.append(f"{attr}={ref}")

    assert missing == []


def test_ipd_app_entrypoints_are_available():
    html = (ROOT / "index.html").read_text(encoding="utf-8")

    assert "e156-submission/assets/ipd-meta-pro.html" in html
    assert (ROOT / "ipd-meta-pro.html").exists()
    assert (ROOT / "e156-submission" / "assets" / "ipd-meta-pro.html").exists()
