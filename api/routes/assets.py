"""Asset management routes."""

from pathlib import Path
from typing import TYPE_CHECKING

from fastapi import APIRouter, File, HTTPException, UploadFile

from lib.exceptions import FileProcessingError

from ..schemas import ListAssetsResponse, UploadResponse

if TYPE_CHECKING:
    from ..config import Settings


def _build_asset_router(settings: "Settings") -> APIRouter:
    """Build the asset router."""

    router = APIRouter(tags=["assets"])

    def _save_upload(target_root: Path, file: UploadFile) -> UploadResponse:
        """Save an uploaded file to the target directory."""
        target_root.mkdir(parents=True, exist_ok=True)
        target_path = target_root / file.filename
        content = file.file.read()
        target_path.write_bytes(content)
        return UploadResponse(name=file.filename, path=str(target_path))

    @router.post("/classifiers", response_model=UploadResponse)
    def upload_classifier(file: UploadFile = File(...)) -> UploadResponse:
        """Upload a classifier file."""
        try:
            return _save_upload(settings.classifier_root, file)
        except FileProcessingError as exc:
            raise HTTPException(
                status_code=500,
                detail=f"File upload error: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Unexpected error: {str(exc)}",
            ) from exc

    @router.get("/classifiers", response_model=ListAssetsResponse)
    def list_classifiers() -> ListAssetsResponse:
        """List all available classifiers."""
        items = [p.name for p in settings.classifier_root.glob("*") if p.is_file()]
        return ListAssetsResponse(items=items)

    @router.post("/adapters", response_model=UploadResponse)
    def upload_adapter(file: UploadFile = File(...)) -> UploadResponse:
        """Upload an adapter file."""
        try:
            return _save_upload(settings.adapter_root, file)
        except FileProcessingError as exc:
            raise HTTPException(
                status_code=500,
                detail=f"File upload error: {exc.message}",
            ) from exc
        except Exception as exc:
            raise HTTPException(
                status_code=500,
                detail=f"Unexpected error: {str(exc)}",
            ) from exc

    @router.get("/adapters", response_model=ListAssetsResponse)
    def list_adapters() -> ListAssetsResponse:
        """List all available adapters."""
        items = [p.name for p in settings.adapter_root.glob("*") if p.is_file()]
        return ListAssetsResponse(items=items)

    @router.post("/delete/classifiers/{name}", response_model=ListAssetsResponse)
    def delete_classifier(name: str) -> ListAssetsResponse:
        """Delete a classifier file."""
        target = settings.classifier_root / name
        if target.exists():
            target.unlink()
        items = [p.name for p in settings.classifier_root.glob("*") if p.is_file()]
        return ListAssetsResponse(items=items)

    @router.post("/delete/adapters/{name}", response_model=ListAssetsResponse)
    def delete_adapter(name: str) -> ListAssetsResponse:
        """Delete an adapter file."""
        target = settings.adapter_root / name
        if target.exists():
            target.unlink()
        items = [p.name for p in settings.adapter_root.glob("*") if p.is_file()]
        return ListAssetsResponse(items=items)

    return router


def register_asset_routes(app, settings: "Settings", *, prefix: str = "", include_in_schema: bool = True) -> None:
    """Register all asset-related routes with the FastAPI app."""
    app.include_router(
        _build_asset_router(settings),
        prefix=prefix,
        include_in_schema=include_in_schema,
    )
