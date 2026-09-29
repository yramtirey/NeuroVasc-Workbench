import { useEffect, useRef, useState } from "react";

import "@kitware/vtk.js/Rendering/Profiles/Geometry";
import vtkActor from "@kitware/vtk.js/Rendering/Core/Actor";
import vtkMapper from "@kitware/vtk.js/Rendering/Core/Mapper";
import vtkGenericRenderWindow from "@kitware/vtk.js/Rendering/Misc/GenericRenderWindow";
import vtkXMLPolyDataReader from "@kitware/vtk.js/IO/XML/XMLPolyDataReader";
import vtkSphereSource from "@kitware/vtk.js/Filters/Sources/SphereSource";
import type vtkPolyData from "@kitware/vtk.js/Common/DataModel/PolyData";

import { getGeometryUrl } from "./api";
import type { CaliberSample } from "./types";

interface VesselViewerProps {
  caseId: string;
  selectedLabel: number;
  labels: number[];
  markerPoint?: CaliberSample | null;
}

interface VesselSurface {
  actor: vtkActor;
  mapper: vtkMapper;
  reader: vtkXMLPolyDataReader;
  polyData: vtkPolyData;
}

interface ViewerScene {
  genericRenderWindow: vtkGenericRenderWindow;
  surfaces: Map<number, VesselSurface>;
  marker: { actor: vtkActor; mapper: vtkMapper; source: vtkSphereSource };
}

function styleActor(actor: vtkActor, selected: boolean) {
  const property = actor.getProperty();
  if (selected) property.setColor(0.10, 0.40, 0.43);
  else property.setColor(0.66, 0.72, 0.73);
  property.setOpacity(selected ? 1 : 0.22);
}

function resetSceneCamera(scene: ViewerScene) {
  const renderer = scene.genericRenderWindow.getRenderer();
  const camera = renderer.getActiveCamera();
  camera.setPosition(0, 0, 1);
  camera.setFocalPoint(0, 0, 0);
  camera.setViewUp(0, 1, 0);
  renderer.resetCamera();
  renderer.resetCameraClippingRange();
  scene.genericRenderWindow.getRenderWindow().render();
}

function VesselViewer({ caseId, selectedLabel, labels, markerPoint = null }: VesselViewerProps) {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const vtkRef = useRef<ViewerScene | null>(null);
  const selectedLabelRef = useRef(selectedLabel);
  const [loadedCount, setLoadedCount] = useState(0);
  const [loadError, setLoadError] = useState<string | null>(null);

  // Keep asynchronous mesh arrivals in sync with the latest committed selection.
  useEffect(() => {
    selectedLabelRef.current = selectedLabel;
    const scene = vtkRef.current;
    if (!scene) return;
    scene.surfaces.forEach(({ actor }, label) => styleActor(actor, label === selectedLabel));
    scene.genericRenderWindow.getRenderWindow().render();
  }, [selectedLabel]);

  // Depend on label values, so an equivalent new array cannot rebuild the scene.
  const labelsKey = JSON.stringify([...new Set(labels)].sort((a, b) => a - b));

  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const controller = new AbortController();
    const sceneLabels: number[] = JSON.parse(labelsKey);
    setLoadedCount(0);
    setLoadError(null);

    const genericRenderWindow = vtkGenericRenderWindow.newInstance({
      background: [0.965, 0.98, 0.98],
      listenWindowResize: false,
    });
    genericRenderWindow.setContainer(container);
    genericRenderWindow.resize();

    const renderer = genericRenderWindow.getRenderer();
    const renderWindow = genericRenderWindow.getRenderWindow();
    const markerSource = vtkSphereSource.newInstance({
      radius: 0.7, thetaResolution: 24, phiResolution: 24,
    });
    const markerMapper = vtkMapper.newInstance({ scalarVisibility: false });
    markerMapper.setInputConnection(markerSource.getOutputPort());
    const markerActor = vtkActor.newInstance();
    markerActor.setMapper(markerMapper);
    markerActor.setVisibility(false);
    markerActor.setPickable(false);
    markerActor.setUseBounds(false); // Reset view should frame anatomy, not the marker.
    markerActor.getProperty().setColor(0.72, 0.30, 0.22);
    markerActor.getProperty().setAmbient(0.45);
    markerActor.getProperty().setDiffuse(0.55);
    renderer.addActor(markerActor);
    const scene: ViewerScene = {
      genericRenderWindow, surfaces: new Map(),
      marker: { actor: markerActor, mapper: markerMapper, source: markerSource },
    };
    vtkRef.current = scene;

    const resizeObserver = new ResizeObserver(() => {
      if (!controller.signal.aborted) genericRenderWindow.resize();
    });
    resizeObserver.observe(container);

    async function loadGeometry() {
      try {
        for (const label of sceneLabels) {
          const response = await fetch(getGeometryUrl(caseId, label), {
            signal: controller.signal,
          });
          if (!response.ok) {
            throw new Error(`Unable to load vessel label ${label} (HTTP ${response.status})`);
          }
          const buffer = await response.arrayBuffer();
          // Aborting fetch alone is insufficient if its response already resolved.
          if (controller.signal.aborted) return;

          const reader = vtkXMLPolyDataReader.newInstance();
          let polyData: vtkPolyData;
          try {
            reader.parseAsArrayBuffer(buffer);
            polyData = reader.getOutputData(0);
            if (!polyData || polyData.getNumberOfPoints() === 0 || polyData.getNumberOfPolys() === 0) {
              throw new Error(`Vessel label ${label} has no surface geometry`);
            }
          } catch (error) {
            reader.getOutputData(0)?.delete();
            reader.delete();
            throw error;
          }

          const mapper = vtkMapper.newInstance({ scalarVisibility: false });
          mapper.setInputData(polyData);
          const actor = vtkActor.newInstance();
          actor.setMapper(mapper);
          const property = actor.getProperty();
          property.setInterpolationToPhong();
          property.setAmbient(0.15);
          property.setDiffuse(0.8);
          property.setSpecular(0.25);
          property.setSpecularPower(15);
          styleActor(actor, label === selectedLabelRef.current);

          scene.surfaces.set(label, { actor, mapper, reader, polyData });
          renderer.addActor(actor);
          setLoadedCount(scene.surfaces.size);
        }
        if (!controller.signal.aborted) resetSceneCamera(scene);
      } catch (error) {
        if (controller.signal.aborted) return;
        setLoadError(error instanceof Error ? error.message : "Unable to load vascular geometry");
      }
    }
    void loadGeometry();

    return () => {
      // Stop asynchronous work and interaction before releasing the OpenGL view.
      controller.abort();
      resizeObserver.disconnect();
      if (vtkRef.current === scene) vtkRef.current = null;
      const interactor = genericRenderWindow.getInteractor();
      const interactorStyle = interactor.getInteractorStyle();
      // vtk.js 37 supports null to detach, but its declaration requires HTMLElement.
      // @ts-expect-error The runtime explicitly handles a null container.
      genericRenderWindow.setContainer(null);
      interactor.delete();
      interactorStyle.delete();
      scene.surfaces.forEach(({ actor, mapper, reader, polyData }) => {
        renderer.removeActor(actor);
        actor.delete();
        mapper.delete();
        reader.delete();
        polyData.delete();
      });
      scene.surfaces.clear();
      renderer.removeActor(markerActor);
      markerActor.getProperty().delete();
      markerActor.delete();
      markerMapper.delete();
      markerSource.getOutputData()?.delete();
      markerSource.delete();
      genericRenderWindow.delete();
      renderer.delete();
      renderWindow.delete();
    };
  }, [caseId, labelsKey]);

  // Update only the existing marker; never reload meshes or move the camera.
  // Scene identity dependencies also apply a supplied point after scene creation.
  useEffect(() => {
    const scene = vtkRef.current;
    if (!scene) return;
    const { actor, source } = scene.marker;
    const valid = markerPoint && [
      markerPoint.x_mm, markerPoint.y_mm, markerPoint.z_mm,
      markerPoint.raw_diameter_mm,
    ].every(Number.isFinite);
    if (valid) {
      source.setCenter(markerPoint.x_mm, markerPoint.y_mm, markerPoint.z_mm);
      // A centerline point lies inside the opaque vessel. Extend just past its
      // local radius so the marker is visible, without changing vessel opacity.
      source.setRadius(Math.max(0.7, markerPoint.raw_diameter_mm / 2 + 0.75));
    }
    actor.setVisibility(!!valid);
    scene.genericRenderWindow.getRenderWindow().render();
  }, [markerPoint, caseId, labelsKey]);

  const totalCount = new Set(labels).size;
  const isLoading = loadedCount < totalCount && !loadError;

  return (
    <div className="vtk-wrapper">
      <div ref={containerRef} className="vtk-container" />
      <div className="viewer-toolbar">
        <button
          type="button"
          className="viewer-tool-button"
          onClick={() => { if (vtkRef.current) resetSceneCamera(vtkRef.current); }}
          disabled={isLoading || !!loadError || loadedCount === 0}
        >
          Reset view
        </button>
      </div>
      {isLoading && (
        <div className="viewer-loading" role="status">
          <div className="loading-spinner" />
          <strong>Loading vascular geometry</strong>
          <span>{loadedCount} / {totalCount} vessels</span>
        </div>
      )}
      {loadError && (
        <div className="viewer-error" role="alert">
          <strong>Unable to load 3D anatomy</strong>
          <span>{loadError}</span>
        </div>
      )}
    </div>
  );
}

export default VesselViewer;
