"""
Custom exception hierarchy for MetaWorks pipeline.

This module defines domain-specific exceptions that provide better error
context and enable more precise error handling throughout the codebase.
"""

from typing import Any, Optional


class MetaWorksError(Exception):
    """Base exception for all MetaWorks-related errors."""

    def __init__(self, message: str, suggestion: Optional[str] = None, **context: Any):
        """
        Initialize MetaWorksError.

        Args:
            message: Error message describing what went wrong
            suggestion: Optional suggestion for how to fix the issue
        """
        self.message = message
        self.suggestion = suggestion
        self.context = context
        super().__init__(self.message)

    def __str__(self):
        parts = [self.message]
        details = self.context.get("details")
        if details:
            parts.append(f"Details: {details}")
        if self.suggestion:
            parts.append(f"Suggestion: {self.suggestion}")
        return "\n".join(parts)


class ConfigurationError(MetaWorksError):
    """Raised when there are issues with configuration."""

    def __init__(
        self,
        message: str,
        config_key: Optional[str] = None,
        filepath: Optional[str] = None,
        suggestion: Optional[str] = None,
        details: Optional[str] = None,
    ):
        """
        Initialize ConfigurationError.

        Args:
            message: Error message
            config_key: Optional key that caused the error
            suggestion: Optional fix suggestion
        """
        self.config_key = config_key
        self.filepath = filepath
        super().__init__(
            message,
            suggestion,
            config_key=config_key,
            filepath=filepath,
            details=details,
        )


class ValidationError(MetaWorksError):
    """Raised when input validation fails."""

    def __init__(
        self,
        message: str,
        field: Optional[str] = None,
        value: Optional[Any] = None,
        suggestion: Optional[str] = None,
        details: Optional[str] = None,
    ):
        """
        Initialize ValidationError.

        Args:
            message: Error message
            field: Optional field name that failed validation
            value: Optional value that failed validation
        """
        self.field = field
        self.value = value
        if suggestion is None and field:
            suggestion = f"Check the '{field}' field"
        super().__init__(message, suggestion, field=field, value=value, details=details)


class FileProcessingError(MetaWorksError):
    """Raised when file processing operations fail."""

    def __init__(
        self,
        message: str,
        filepath: Optional[str] = None,
        suggestion: Optional[str] = None,
        details: Optional[str] = None,
    ):
        """
        Initialize FileProcessingError.

        Args:
            message: Error message
            filepath: Optional path to the problematic file
            suggestion: Optional fix suggestion
        """
        self.filepath = filepath
        super().__init__(message, suggestion, filepath=filepath, details=details)


class PipelineExecutionError(MetaWorksError):
    """Raised when pipeline execution fails."""

    def __init__(
        self,
        message: str,
        stage: Optional[str] = None,
        exit_code: Optional[int] = None,
        suggestion: Optional[str] = None,
        details: Optional[str] = None,
    ):
        """
        Initialize PipelineExecutionError.

        Args:
            message: Error message
            stage: Optional pipeline stage that failed
            exit_code: Optional exit code from the failed process
            suggestion: Optional fix suggestion
        """
        self.stage = stage
        self.exit_code = exit_code
        super().__init__(message, suggestion, stage=stage, exit_code=exit_code, details=details)


class MetaWorksRuntimeError(MetaWorksError):
    """Raised when runtime-specific operations fail."""

    def __init__(
        self,
        message: str,
        runtime_type: Optional[str] = None,
        suggestion: Optional[str] = None,
        details: Optional[str] = None,
    ):
        self.runtime_type = runtime_type
        super().__init__(message, suggestion, runtime_type=runtime_type, details=details)


class StateManagementError(MetaWorksError):
    """Raised when state management operations fail."""

    def __init__(
        self,
        message: str,
        run_id: Optional[str] = None,
        suggestion: Optional[str] = None,
        details: Optional[str] = None,
    ):
        """
        Initialize StateManagementError.

        Args:
            message: Error message
            run_id: Optional run ID that caused the error
            suggestion: Optional fix suggestion
        """
        self.run_id = run_id
        super().__init__(message, suggestion, run_id=run_id, details=details)


class DependencyError(MetaWorksError):
    """Raised when required dependencies are missing or incompatible."""

    def __init__(
        self,
        message: str,
        dependency: Optional[str] = None,
        suggestion: Optional[str] = None,
        details: Optional[str] = None,
    ):
        """
        Initialize DependencyError.

        Args:
            message: Error message
            dependency: Optional missing dependency name
            suggestion: Optional fix suggestion
        """
        self.dependency = dependency
        if not suggestion and dependency:
            suggestion = f"Install {dependency}"
        super().__init__(message, suggestion, dependency=dependency, details=details)


class DataError(MetaWorksError):
    """Raised when data processing issues occur."""

    def __init__(
        self,
        message: str,
        data_type: Optional[str] = None,
        suggestion: Optional[str] = None,
        details: Optional[str] = None,
    ):
        """
        Initialize DataError.

        Args:
            message: Error message
            data_type: Optional type of data that caused the error
            suggestion: Optional fix suggestion
        """
        self.data_type = data_type
        super().__init__(message, suggestion, data_type=data_type, details=details)
