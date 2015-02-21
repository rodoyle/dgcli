# DeskGen Platform Command Line Interface

## Design Library

A library is a single task, that should spawn one sub-task per target. These targets should be designed
in parallel. After each sub-task completes the resulting library should be combined and saved to a single
notebook.

The generally calling conventions should follow JSON RPC and the output shall be a JSON document.
Complex library designs should be specified by nested JSON RPC calls to the
relevant subtasks.

Automatically added elements are indicated with AUTO

Example Library Input
```javascript
{
  "jsonrpc": '2.0', //AUTO
  "method": 'design_library', // Always use a verb
  "id": 1, //AUTO
  "params": {
    "genome": {'version': "GRCh38.p1"}, // Always use an object as an identifier
    "nuclease": {"name": "wtCas9"},
    "defaults": { 
      // parameters object, the default design rules for each target
      "method": 'walk_gene', // or similar function
      "params": {
        // Default Parameters Object
      }
    },
    "targets": [
       // an array of targets, each element should be a JSON RPC request object
      {
        "scoring_function": { a different scoring function},
        "gene": {gene dictionary},
        "filters": [ list of filters],
      },
    ],
    "callbacks": [
      // An array of post-parallel execution tasks to call after all
      // targets are run through the system, such as save to notebook,
      // Might not keep this
      {
        "method": 'save_to_notebook',
      },
      {
        "methods": 'save_to_file',
      }
    ],
    "user": {
      // User object containing the email, Auth token, etc for the user
    },
    
    "name": "MouseKnockout Library, // A name for the library
    "description": "A test library for knocking out mouse genes",
    "date_created": date stamp 
    "format": {
      // Formatting speficications
    },
    "dry_run": true, //Don't save anything
    "async": false, // Return a task ID or wait for the output stream,
  }
}    
``` 

Library Result Object

{
  "id": 1,
  "result": {
    "name": 
    "guides": [
        // An array of results for each target in the library
        // Each element either a guide, pair, or vector
        // Could be grouped by target, but probably cleaner to have each element
        // Be a distinct molecular entity. 
    ], 
  }
}