from pathlib import Path
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import subprocess
from threading import Thread

import pytest
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait


ROOT = Path(__file__).resolve().parents[1]


class QuietHTTPServer(ThreadingHTTPServer):
    allow_reuse_address = False
    daemon_threads = True


class QuietHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT), **kwargs)

    def log_message(self, format, *args):
        return


@pytest.fixture(scope="module")
def base_url():
    server = QuietHTTPServer(("127.0.0.1", 0), QuietHandler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_port}"
    finally:
        server.shutdown()
        thread.join(timeout=5)
        server.server_close()


@pytest.fixture()
def driver(base_url):
    options = Options()
    options.add_argument("--headless=new")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    browser = webdriver.Chrome(options=options)
    try:
        yield browser
    finally:
        browser.quit()


def test_index_page_loads(driver, base_url):
    driver.get(f"{base_url}/index.html")
    WebDriverWait(driver, 10).until(lambda d: d.find_element(By.TAG_NAME, "h1").text.strip() != "")

    assert "IPD Meta-Analysis Pro" in driver.title
    assert driver.find_element(By.TAG_NAME, "h1").text.startswith("IPD Meta-Analysis Pro")
    links = driver.find_elements(By.CSS_SELECTOR, ".links a.card")
    assert len(links) >= 4
    figures = driver.find_elements(By.CSS_SELECTOR, ".figs img")
    assert len(figures) >= 6


def test_single_file_app_loads(driver, base_url):
    driver.get(f"{base_url}/ipd-meta-pro.html")
    WebDriverWait(driver, 10).until(lambda d: d.find_element(By.CLASS_NAME, "logo-text").text.strip() != "")

    assert driver.title.startswith("IPD Meta-Analysis Pro")
    assert driver.find_element(By.CLASS_NAME, "logo-text").text == "IPD Meta-Analysis Pro"
    assert driver.find_element(By.ID, "sessionLabel").text.strip()
    assert len(driver.find_elements(By.CSS_SELECTOR, ".nav-tab")) >= 4


def test_core_node_suite_passes():
    completed = subprocess.run(
        ["node", "dev/tests/core_stats_meta_test.js"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Core stats/meta module test passed" in completed.stdout
