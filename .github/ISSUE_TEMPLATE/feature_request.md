---
name: Feature request
about: Suggest a new marker, module, or workflow capability
title: "[feature] "
labels: enhancement
body:
  - type: textarea
    id: problem
    attributes:
      label: Problem or motivation
      description: What are you trying to do that MetaWorks does not support today?
    validations:
      required: true
  - type: textarea
    id: proposal
    attributes:
      label: Proposed solution
      description: What should change? A new module, preset, parameter, or UI capability?
    validations:
      required: true
  - type: textarea
    id: alternatives
    attributes:
      label: Alternatives considered
      description: Existing tools or current workarounds you have tried
