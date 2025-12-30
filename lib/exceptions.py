"""
Custom exception hierarchy for MetaWorks pipeline.

This module defines domain-specific exceptions that provide better error
context and enable more precise error handling throughout the codebase.
"""


class MetaWorksError(Exception):
    """Base exception for all MetaWorks-related errors."""

    def __init__(self, message: str, suggestion: str = None):
        """
        Initialize MetaWorksError.

        Args:
            message: Error message describing what went wrong
            suggestion: Optional suggestion for how to fix the issue
        """
        self.message = message
        self.suggestion = suggestion
        super().__init__(self.message)

    def __str__(self):
        if self.suggestion:
            return f"{self.message}\nSuggestion: {self.suggestion}"
        return self.message


class ConfigurationError(MetaWorksError):
    """Raised when there are issues with configuration."""

    def __init__(self, message: str, config_key: str = None, suggestion: str = None):
        """
        Initialize ConfigurationError.

        Args:
            message: Error message
            config_key: Optional key that caused the error
            suggestion: Optional fix suggestion
        """
        self.config_key = config_key
        super().__init__(message, suggestion)


class ValidationError(MetaWorksError):
    """Raised when input validation fails."""

    def __init__(self, message: str, field: str = None, value: any = None):
        """
        Initialize ValidationError.

        Args:
            message: Error message
            field: Optional field name that failed validation
            value: Optional value that failed validation
        """
        self.field = field
        self.value = value
        suggestion = f"Check the '{field}' field" if field else None
        super().__init__(message, suggestion)


class FileProcessingError(MetaWorksError):
    """Raised when file processing operations fail."""

    def __init__(self, message: str, filepath: str = None, suggestion: str = None):
        """
        Initialize FileProcessingError.

        Args:
            message: Error message
            filepath: Optional path to the problematic file
            suggestion: Optional fix suggestion
        """
        self.filepath = filepath
        super().__init__(message, suggestion)


class PipelineExecutionError(MetaWorksError):
    """Raised when pipeline execution fails."""

    def __init__(
        self, message: str, stage: str = None, exit_code: int = None, suggestion: str = None
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
        super().__init__(message, suggestion)


class RuntimeError(MetaWorksError):
    """Raised when runtime-specific operations fail."""

    def __init__(self, message: str, runtime_type: str = None, suggestion: str = None):
        """
        Initialize RuntimeError.

        Args:
            message: Error message
            runtime_type: Optional runtime type (conda, docker, singularity)
            suggestion: Optional fix suggestion
        """
        self.runtime_type = runtime_type
        super().__init__(message, suggestion)


class StateManagementError(MetaWorksError):
    """Raised when state management operations fail."""

    def __init__(self, message: str, run_id: str = None, suggestion: str = None):
        """
        Initialize StateManagementError.

        Args:
            message: Error message
            run_id: Optional run ID that caused the error
            suggestion: Optional fix suggestion
        """
        self.run_id = run_id
        super().__init__(message, suggestion)


class DependencyError(MetaWorksError):
    """Raised when required dependencies are missing or incompatible."""

    def __init__(self, message: str, dependency: str = None, suggestion: str = None):
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
        super().__init__(message, suggestion)


class DataError(MetaWorksError):
    """Raised when data processing issues occur."""

    def __init__(self, message: str, data_type: str = None, suggestion: str = None):
        """
        Initialize DataError.

        Args:
            message: Error message
            data_type: Optional type of data that caused the error
            suggestion: Optional fix suggestion
        """
        self.data_type = data_type
        super().__init__(message, suggestion)
