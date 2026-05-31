"""
Rate limiting and retry utilities for API clients.

UniProt has rate limits. This module provides a simple
throttle + exponential backoff wrapper so API clients don't have to
implement it themselves.
"""
from __future__ import annotations

import logging
import time
from typing import Callable, TypeVar

import requests

log = logging.getLogger(__name__)

T = TypeVar("T")


def get_with_retry(
    url: str,
    params: dict | None = None,
    headers: dict | None = None,
    max_retries: int = 5,
    backoff_base: float = 2.0,
    timeout: int = 60,
) -> requests.Response:
    """GET with exponential backoff on 429 / 5xx responses."""
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, params=params, headers=headers, timeout=timeout)
            if resp.status_code == 200:
                return resp
            if resp.status_code == 429:
                wait = backoff_base ** attempt
                log.warning("Rate limited (429) on %s — waiting %.1fs", url, wait)
                time.sleep(wait)
                continue
            if resp.status_code >= 500:
                wait = backoff_base ** attempt
                log.warning("Server error %d on %s — waiting %.1fs", resp.status_code, url, wait)
                time.sleep(wait)
                continue
            # 4xx that isn't 429 — no point retrying
            resp.raise_for_status()
        except requests.exceptions.Timeout:
            wait = backoff_base ** attempt
            log.warning("Timeout on %s — waiting %.1fs", url, wait)
            time.sleep(wait)
        except requests.exceptions.ConnectionError as e:
            wait = backoff_base ** attempt
            log.warning("Connection error on %s: %s — waiting %.1fs", url, e, wait)
            time.sleep(wait)

    raise RuntimeError(f"Failed to GET {url} after {max_retries} attempts")


def post_with_retry(
    url: str,
    data: dict | None = None,
    json: dict | None = None,
    headers: dict | None = None,
    max_retries: int = 5,
    backoff_base: float = 2.0,
    timeout: int = 120,
) -> requests.Response:
    """POST with exponential backoff on 429 / 5xx responses."""
    for attempt in range(max_retries):
        try:
            resp = requests.post(url, data=data, json=json, headers=headers, timeout=timeout)
            if resp.status_code == 200:
                return resp
            if resp.status_code == 429:
                wait = backoff_base ** attempt
                log.warning("Rate limited (429) — waiting %.1fs", wait)
                time.sleep(wait)
                continue
            if resp.status_code >= 500:
                wait = backoff_base ** attempt
                log.warning("Server error %d — waiting %.1fs", resp.status_code, wait)
                time.sleep(wait)
                continue
            resp.raise_for_status()
        except requests.exceptions.Timeout:
            wait = backoff_base ** attempt
            log.warning("Timeout — waiting %.1fs", wait)
            time.sleep(wait)
        except requests.exceptions.ConnectionError as e:
            wait = backoff_base ** attempt
            log.warning("Connection error: %s — waiting %.1fs", e, wait)
            time.sleep(wait)

    raise RuntimeError(f"Failed to POST {url} after {max_retries} attempts")
