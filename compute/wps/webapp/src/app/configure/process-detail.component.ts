import { Component, OnInit, Input, Output, EventEmitter } from '@angular/core';

import { Dataset } from './dataset';
import { Process } from './process';
import { Variable, VariableCollection, File } from './configure.service';
import { RegridModel } from './regrid.component';
import { Parameter } from './parameter.component';
import { Axis } from './domain.component';
import { Domain } from './domain.component';

import { AuthService } from '../core/auth.service';
import { ConfigureService } from './configure.service';
import { ConfigService } from '../core/config.service';
import { WPSService } from '../core/wps.service';
import { NotificationService } from '../core/notification.service';

import { Subject } from 'rxjs/Subject';

declare var $: any;

class FileModalState {
  constructor(
    public n: number,
    public start: number = 0,
    public stop: number = 20,
    public step: number = 20,
  ) { }

  canPrevious() {
    return (this.start-this.step) >= 0;
  }

  canNext() {
    return (this.stop+this.step) <= this.n;
  }

  previous() {
    if (!this.canPrevious()) { return; }

    this.start = Math.max(this.start-this.step, 0);

    this.stop -= this.step;
  }

  next() {
    if (!this.canNext()) { return; }

    this.start += this.step;

    this.stop = Math.min(this.stop+this.step, this.n);
  }
}

@Component({
  selector: 'process-detail',
  styles: [`
  .file-table-container {
    height: 70vh;
    overflow: auto;
  }

  .file-table {
    table-layout: fixed;
  }

  .file-url-column {
    width: 90%;
    word-wrap: break-word;
  }
  `],
  templateUrl: './process-detail.component.html',
  providers: [ConfigureService],
})
export class ProcessDetailComponent implements OnInit {
  @Input() params: any;
  @Output() modified = new EventEmitter<Domain>();

  _datasets: Dataset[];
  operations: Process[];
  variables: Variable[];

  selectedDataset: Dataset;
  selectedOperation: Process;
  selectedVariable: Variable;

  files: File[] = [];

  regridModel = new RegridModel('ESMF', 'Linear', 'None', -90, 180, 0.5, 0, 360, 0.5);
  parameters: Parameter[] = [];
  axes: Axis[] = [];

  fileState: FileModalState;

  filterPattern = '';
  _filterPattern = new Subject<string>();

  constructor(
    public authService: AuthService,
    public configService: ConfigService,
    public configureService: ConfigureService,
    public notificationService: NotificationService,
    public wpsService: WPSService,
  ) { }

  ngOnInit() {
    this.wpsService.getCapabilities(this.configService.wpsPath).then((processes: Process[]) => {
      this.operations = processes.sort((a: Process, b: Process) => {
        if (a.identifier < b.identifier) { return -1; }
        if (a.identifier > b.identifier) { return 1; }

        return 0;
      });

      this.changeOperation(this.operations[0]);
    }).catch((error: any) => {
      this.notificationService.error(error); 
    });

    this._filterPattern.
      debounceTime(100).
      distinctUntilChanged().
      subscribe((value: string) => {
        this.filterPattern = value;
      });
  }

  @Input()
  set datasets(values: string[]) {
    this._datasets = values.map((id: string) => {
      return new Dataset(id);
    });

    if (this._datasets.length > 0) {
      this.changeDataset(this._datasets[0]);
    }
  }

  downloadScript() {
    let process = this.selectedOperation;

    process.inputs = this.files.filter((item: File) => { return item.included; });

    if (this.selectedVariable != null) {
      process.variable = this.selectedVariable.id;
    }

    process.domain = new Domain('', this.axes);

    process.regrid = this.regridModel;

    process.parameters = this.parameters;

    this.configureService.downloadScript(process).then((data: any) => {
      let url = URL.createObjectURL(new Blob([data.text]));

      let a = document.createElement('a');

      a.href = url;
      a.target = '_blank';
      a.download = data.filename;

      a.click();
    }).catch((error: any) => {
      this.notificationService.error(error); 
    });
  }

  getSelectedFiles() {
    if (this.files == null) { throw 'No files have been selected'; }

    return this.files.filter((item: File) => { return item.included; });
  }

  getSelectedVariableName() {
    if (this.selectedVariable == null) { throw 'No Variable has been selected'; }

    return this.selectedVariable.id;
  }

  execute() {
    let process = this.selectedOperation;

    if (process.description.metadata.inputs > 0) {
      try {
        process.inputs = this.getSelectedFiles();

        process.variable = this.getSelectedVariableName();

        process.domain = new Domain('', this.axes);

        process.regrid = this.regridModel;
      } catch(e) {
        this.notificationService.error(e);

        return;
      }
    }

    process.parameters = this.parameters;

    this.wpsService.execute(
      this.configService.wpsPath,
      this.authService.user.api_key,
      process
    ).then((data: any) => {
      this.notificationService.message(data);      
    }).catch((error: any) => {
      this.notificationService.error(error);
    });
  }

  showFiles() {
    $('#filesModal').on('hidden.bs.modal', (event: any) => {
      $('#filesModal').off('hidden.bs.modal');

      this.updateVariableFiles();
    });
  }

  changeVariable(variable: Variable) {
    this.selectedVariable = variable;

    let result = this.selectedDataset.variableCollection;

    this.files = this.selectedVariable.files.map((key: number, index: number) => {
      return new File(`${result.files[key]}`, key, (index == 0 ? true : false));
    });

    this.fileState = new FileModalState(this.files.length);

    this.files.sort((a: File, b: File) => {
      if (a.url < b.url) { return -1; }
      if (a.url > b.url) { return 1; }

      return 0;
    });

    this.updateVariableFiles();
  }

  updateVariableFiles() {
    let selectedFiles = this.files.filter((item: File) => {
      return item.included;
    });

    if (selectedFiles.length == 0) {
      this.notificationService.warn('You have not files selected');

      this.axes = [];

      return;
    }
    
    let selectedIndexes = selectedFiles.map((item: File) => {
      return item.index; 
    });

    this.configureService.searchVariable(
      this.selectedVariable.id,
      this.selectedDataset.masterID,
      selectedIndexes,
      this.params,
    ).then((data: any[]) => {
      for (let x = 0; x < data.length; x++) {
        let file = this.files.find((item: File) => {
          return item.url == data[x].url;
        });

        file.temporal = data[x].temporal;

        file.spatial = data[x].spatial;
      }

      this.updateDomain(selectedFiles);
    }).catch((error: any) => {
      this.notificationService.error(error); 
    });
  }

  updateDomain(selectedFiles: File[]) {
    let spatial = null;
    let temporal = null;

    for (let x = 0; x < selectedFiles.length; x++) {
      let file = selectedFiles[x];

      if (spatial == null) {
        spatial = [];

        for (let y = 0; y < file.spatial.length; y++) {
          let axis = file.spatial[y];

          let newAxis = new Axis(
            axis.id,
            axis.units,
            axis.start,
            axis.stop,
            axis.length,
            'spatial',
          );

          newAxis.data = axis;

          spatial.push(newAxis);
        }
      }

      if (temporal == null) {
        temporal = new Axis(
          file.temporal.id,
          file.temporal.units,
          file.temporal.start,
          file.temporal.stop,
          file.temporal.length,
          'temporal',
        );

        temporal.data = selectedFiles;
      } else {
        if (temporal._start > file.temporal.start) {
          temporal.start = temporal._start = file.temporal.start;
        }

        if (temporal._stop < file.temporal.stop) {
          temporal.stop = temporal._stop = file.temporal.stop;
        }
      }
    }

    this.axes.length = 0;

    this.axes.push(temporal);
    this.axes.push(...spatial);
  }

  changeDataset(dataset: Dataset) {
    this.selectedDataset = dataset;

    this.configureService.searchESGF(
      this.selectedDataset.masterID, 
      this.params
    ).then((result: VariableCollection) => {
      result.variables.sort((a: Variable, b: Variable) => {
        if (a.id < b.id) { return -1; }
        if (a.id > b.id) { return 1; }

        return 0;
      });

      this.selectedDataset.variableCollection = result;

      this.changeVariable(result.variables[0]);
    }).catch((error: any) => {
      this.notificationService.error(error); 
    });
  }

  changeOperation(operation: Process) {
    this.selectedOperation = operation;

    if (this.selectedOperation.description == null) {
      this.wpsService.describeProcess(
        this.configService.wpsPath, 
        this.selectedOperation.identifier
      ).then((data: any) => {
        this.selectedOperation.description = data;
      });
    }
  }

  updateFilesIncluded(included: boolean) {
    this.files.forEach((file: File) => {
      file.included = included;
    });
  }

  addParameter() {
    this.parameters.push(new Parameter());
  }

  removeParameter(param: Parameter) {
    this.parameters = this.parameters.filter((item: Parameter) => {
      return item.uid != param.uid;
    });
  }
}
