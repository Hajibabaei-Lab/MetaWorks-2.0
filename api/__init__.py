"""
API package exposing MetaWorks workflows through FastAPI.

This module intentionally keeps initialization light so uvicorn can import
`api.main:app` without executing any heavy logic up front.
"""

__all__ = ["__version__"]

__version__ = "0.1.0"
