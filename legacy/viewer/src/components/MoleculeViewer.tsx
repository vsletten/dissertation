import { useRef, useState, useMemo, useCallback } from 'react';
import { Canvas, useThree } from '@react-three/fiber';
import { OrbitControls, PerspectiveCamera, Line, Text, Billboard } from '@react-three/drei';
import * as THREE from 'three';
import { getElementColor, getCovalentRadius, getVdwRadius } from '../utils/elementColors';

interface Atom {
  id: number;
  element: string;
  position: [number, number, number];
  label: string;
  charge: number;
}

interface Bond {
  from: number;
  to: number;
}

interface MoleculeViewerProps {
  atoms: Atom[];
  bonds: Bond[];
}

// Display options
interface DisplayOptions {
  showLabels: boolean;
  displayStyle: 'ball-and-stick' | 'space-filling' | 'wireframe';
  bondRadius: number;
  atomSizeScale: number;
  colorScheme: 'element' | 'charge';
}

// Controls for camera and interaction
const Controls = ({ controlsRef }: { controlsRef: React.RefObject<any> }) => {
  return (
    <OrbitControls 
      ref={controlsRef}
      enableDamping
      dampingFactor={0.2}
      rotateSpeed={0.8}
      zoomSpeed={1.0}
      panSpeed={0.8}
      minDistance={0.1}
      maxDistance={10000}
    />
  );
};

// Individual atom representation
const AtomSphere = ({
  position,
  element,
  covalentRadius,
  vdwRadius,
  charge,
  label,
  options,
}: {
  position: [number, number, number];
  element: string;
  covalentRadius: number;
  vdwRadius: number;
  charge: number;
  label: string;
  options: DisplayOptions;
}) => {
  // Use element-based coloring or charge-based coloring
  const color = useMemo(() => {
    if (options.colorScheme === 'element') {
      return getElementColor(element);
    } else {
      // Charge-based coloring - blue for positive, red for negative
      if (charge > 0) {
        const intensity = Math.min(charge * 2, 1); // Scale the intensity
        return `rgb(${Math.floor(255 * (1 - intensity))}, ${Math.floor(255 * (1 - intensity))}, 255)`;
      } else if (charge < 0) {
        const intensity = Math.min(Math.abs(charge) * 2, 1); // Scale the intensity
        return `rgb(255, ${Math.floor(255 * (1 - intensity))}, ${Math.floor(255 * (1 - intensity))})`;
      }
      return '#AAAAAA'; // Neutral
    }
  }, [element, charge, options.colorScheme]);

  // Adjust radius based on display style
  // Space-filling uses van der Waals radii (true atomic volume)
  // Ball-and-stick and wireframe use covalent radii (bonding size)
  const displayRadius = useMemo(() => {
    if (options.displayStyle === 'space-filling') {
      return vdwRadius * options.atomSizeScale;
    } else if (options.displayStyle === 'ball-and-stick') {
      return covalentRadius * 0.5 * options.atomSizeScale;
    } else {
      return covalentRadius * 0.3 * options.atomSizeScale;
    }
  }, [covalentRadius, vdwRadius, options.displayStyle, options.atomSizeScale]);

  // Text label for the atom using Billboard to always face camera
  const TextLabel = () => {
    if (!options.showLabels) return null;

    const labelText = label || element;

    return (
      <Billboard position={[position[0], position[1] + displayRadius + 0.3, position[2]]}>
        <Text
          fontSize={0.4}
          color="#ffffff"
          anchorX="center"
          anchorY="middle"
          outlineWidth={0.02}
          outlineColor="#000000"
        >
          {labelText}
        </Text>
      </Billboard>
    );
  };

  // Render based on display style
  return (
    <>
      {options.displayStyle !== 'wireframe' ? (
        <mesh position={position}>
          <sphereGeometry args={[displayRadius, 32, 32]} />
          <meshPhongMaterial color={color} />
        </mesh>
      ) : (
        <mesh position={position}>
          <sphereGeometry args={[displayRadius, 16, 16]} />
          <meshBasicMaterial color={color} wireframe />
        </mesh>
      )}
      <TextLabel />
    </>
  );
};

// Bond representation
const Bond = ({
  start,
  end,
  options,
}: {
  start: [number, number, number];
  end: [number, number, number];
  options: DisplayOptions;
}) => {
  // Don't show bonds in space-filling mode
  if (options.displayStyle === 'space-filling') {
    return null;
  }

  const middlePoint = useMemo(() => {
    return new THREE.Vector3(
      (start[0] + end[0]) / 2,
      (start[1] + end[1]) / 2,
      (start[2] + end[2]) / 2
    );
  }, [start, end]);

  // Calculate bond direction and length
  const { length, quaternion } = useMemo(() => {
    const startVec = new THREE.Vector3(...start);
    const endVec = new THREE.Vector3(...end);
    const direction = new THREE.Vector3().subVectors(endVec, startVec);
    const len = direction.length();

    // Cylinder is created along Y-axis, so we need to rotate from Y to the bond direction
    const yAxis = new THREE.Vector3(0, 1, 0);
    const quat = new THREE.Quaternion();
    quat.setFromUnitVectors(yAxis, direction.clone().normalize());

    return { length: len, quaternion: quat };
  }, [start, end]);

  // Determine bond radius based on display style
  const bondRadius = useMemo(() => {
    if (options.displayStyle === 'wireframe') {
      return 0.02 * options.bondRadius;
    }
    return 0.1 * options.bondRadius;
  }, [options.displayStyle, options.bondRadius]);

  // In wireframe mode, use a simple line
  if (options.displayStyle === 'wireframe') {
    return (
      <Line
        points={[
          new THREE.Vector3(...start),
          new THREE.Vector3(...end)
        ]}
        color="#FFFFFF"
        lineWidth={bondRadius * 100}
      />
    );
  }

  // For ball-and-stick, use a cylinder
  return (
    <mesh position={middlePoint} quaternion={quaternion}>
      <cylinderGeometry args={[bondRadius, bondRadius, length, 12]} />
      <meshPhongMaterial color="#CCCCCC" />
    </mesh>
  );
};

// Main Molecule component
const Molecule = ({ 
  atoms, 
  bonds, 
  options, 
  modelCenter 
}: MoleculeViewerProps & { 
  options: DisplayOptions,
  modelCenter: [number, number, number]
}) => {
  // Create a map of atom IDs to array indices for efficient lookup
  const atomMap = useMemo(() => {
    return atoms.reduce((map, atom, index) => {
      map[atom.id] = index;
      return map;
    }, {} as Record<number, number>);
  }, [atoms]);

  // Calculate centered positions for atoms
  const centeredAtoms = useMemo(() => {
    return atoms.map(atom => ({
      ...atom,
      centeredPosition: [
        atom.position[0] - modelCenter[0],
        atom.position[1] - modelCenter[1],
        atom.position[2] - modelCenter[2]
      ] as [number, number, number]
    }));
  }, [atoms, modelCenter]);

  return (
    <group>
      {/* Render all atoms */}
      {centeredAtoms.map((atom) => (
        <AtomSphere
          key={`atom-${atom.id}`}
          position={atom.centeredPosition}
          element={atom.element}
          covalentRadius={getCovalentRadius(atom.element)}
          vdwRadius={getVdwRadius(atom.element)}
          charge={atom.charge}
          label={atom.label}
          options={options}
        />
      ))}
      
      {/* Render all bonds */}
      {bonds.map((bond, index) => {
        const fromAtom = centeredAtoms[atomMap[bond.from]];
        const toAtom = centeredAtoms[atomMap[bond.to]];
        
        if (!fromAtom || !toAtom) {
          return null; // Skip invalid bonds
        }
        
        return (
          <Bond
            key={`bond-${index}`}
            start={fromAtom.centeredPosition}
            end={toAtom.centeredPosition}
            options={options}
          />
        );
      })}
    </group>
  );
};

// Main viewer component
export default function MoleculeViewer({ atoms, bonds }: MoleculeViewerProps): React.ReactElement {
  const controlsRef = useRef<any>(null);
  const [controlsCollapsed, setControlsCollapsed] = useState(false);

  // Display options with defaults
  const [options, setOptions] = useState<DisplayOptions>({
    showLabels: true,
    displayStyle: 'ball-and-stick',
    bondRadius: 1.0,
    atomSizeScale: 1.0,
    colorScheme: 'element',
  });
  
  // Function to reset the camera view
  const resetView = useCallback(() => {
    if (controlsRef.current) {
      controlsRef.current.reset();
    }
  }, [controlsRef]);

  // Calculate model bounds and center
  const { center, boundingBox } = useMemo(() => {
    if (atoms.length === 0) return { center: [0, 0, 0], boundingBox: { min: [0, 0, 0], max: [0, 0, 0] } };
    
    // Find min and max for each axis
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    let minZ = Infinity, maxZ = -Infinity;
    
    atoms.forEach((atom) => {
      minX = Math.min(minX, atom.position[0]);
      maxX = Math.max(maxX, atom.position[0]);
      minY = Math.min(minY, atom.position[1]);
      maxY = Math.max(maxY, atom.position[1]);
      minZ = Math.min(minZ, atom.position[2]);
      maxZ = Math.max(maxZ, atom.position[2]);
    });
    
    // Calculate center from bounds
    const centerX = (minX + maxX) / 2;
    const centerY = (minY + maxY) / 2;
    const centerZ = (minZ + maxZ) / 2;
    
    return {
      center: [centerX, centerY, centerZ],
      boundingBox: {
        min: [minX, minY, minZ],
        max: [maxX, maxY, maxZ]
      }
    };
  }, [atoms]);

  // Calculate a reasonable camera distance based on model size
  const cameraDistance = useMemo(() => {
    if (atoms.length === 0) return 20;
    
    // Calculate model size
    const size = [
      boundingBox.max[0] - boundingBox.min[0],
      boundingBox.max[1] - boundingBox.min[1],
      boundingBox.max[2] - boundingBox.min[2]
    ];
    
    // Use the largest dimension for the camera distance
    const maxSize = Math.max(...size);
    
    // Add padding factor for better viewing
    // Scale factor is reduced for very large models
    const scaleFactor = maxSize > 500 ? 1.2 : 
                       maxSize > 200 ? 1.3 : 
                       maxSize > 100 ? 1.4 : 1.5;
    
    return maxSize * scaleFactor;
  }, [atoms, boundingBox]);

  // Handle option changes
  const handleOptionChange = useCallback((key: keyof DisplayOptions, value: any) => {
    setOptions(prev => ({ ...prev, [key]: value }));
  }, []);

  return (
    <div className="molecule-viewer">
      {controlsCollapsed ? (
        <button
          className="controls-toggle collapsed"
          onClick={() => setControlsCollapsed(false)}
          title="Show controls"
        >
          <span className="toggle-icon">+</span>
        </button>
      ) : (
        <div className="viewer-controls">
          <button
            className="controls-collapse-btn"
            onClick={() => setControlsCollapsed(true)}
            title="Hide controls"
          >
            -
          </button>
          <div className="control-group">
            <label>
              <input
                type="checkbox"
                checked={options.showLabels}
                onChange={(e) => handleOptionChange('showLabels', e.target.checked)}
              />
              Show Labels
            </label>
          </div>

          <div className="control-group">
            <label>Display Style:</label>
            <select
              value={options.displayStyle}
              onChange={(e) => handleOptionChange('displayStyle', e.target.value as any)}
            >
              <option value="ball-and-stick">Ball and Stick</option>
              <option value="space-filling">Space Filling</option>
              <option value="wireframe">Wireframe</option>
            </select>
          </div>

          <div className="control-group">
            <label>Color By:</label>
            <select
              value={options.colorScheme}
              onChange={(e) => handleOptionChange('colorScheme', e.target.value as any)}
            >
              <option value="element">Element</option>
              <option value="charge">Charge</option>
            </select>
          </div>

          <div className="control-group">
            <label>
              Atom Size: {options.atomSizeScale.toFixed(1)}
            </label>
            <input
              type="range"
              min="0.5"
              max="2.0"
              step="0.1"
              value={options.atomSizeScale}
              onChange={(e) => handleOptionChange('atomSizeScale', parseFloat(e.target.value))}
            />
          </div>

          <div className="control-group">
            <label>
              Bond Size: {options.bondRadius.toFixed(1)}
            </label>
            <input
              type="range"
              min="0.5"
              max="2.0"
              step="0.1"
              value={options.bondRadius}
              onChange={(e) => handleOptionChange('bondRadius', parseFloat(e.target.value))}
            />
          </div>

          <div className="control-group">
            <button
              className="reset-view-button"
              onClick={resetView}
              title="Reset camera view"
            >
              Reset View
            </button>
          </div>
        </div>
      )}
      
      <div className="viewer-container">
        <Canvas dpr={[1, 2]} style={{ background: '#111' }}>
          <PerspectiveCamera
            makeDefault
            position={[0, 0, cameraDistance]}
            fov={45}
            near={0.01}
            far={100000}
          />
          <ambientLight intensity={0.5} />
          <directionalLight position={[10, 10, 10]} intensity={0.7} />
          <directionalLight position={[-10, -10, -10]} intensity={0.5} />
          <Controls controlsRef={controlsRef} />
          {atoms.length > 0 && (
            <Molecule 
              atoms={atoms} 
              bonds={bonds} 
              options={options} 
              modelCenter={center as [number, number, number]} 
            />
          )}
        </Canvas>
      </div>
      
      <style>{`
        .molecule-viewer {
          display: flex;
          flex-direction: column;
          width: 100%;
          height: 100%;
        }

        .viewer-controls {
          position: absolute;
          top: 10px;
          left: 10px;
          display: flex;
          flex-wrap: wrap;
          gap: 10px 16px;
          padding: 10px 14px;
          padding-top: 28px;
          background-color: rgba(10, 20, 30, 0.85);
          backdrop-filter: blur(5px);
          color: #7AFFB2;
          z-index: 10;
          border-radius: 4px;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
          border: 1px solid rgba(7, 219, 250, 0.3);
          box-shadow: 0 0 15px rgba(7, 219, 250, 0.2);
          max-width: 520px;
        }

        .controls-collapse-btn {
          position: absolute;
          top: 4px;
          right: 4px;
          width: 20px;
          height: 20px;
          padding: 0;
          background-color: rgba(7, 219, 250, 0.2);
          color: #07DBFA;
          border: 1px solid rgba(7, 219, 250, 0.4);
          border-radius: 3px;
          font-size: 16px;
          font-weight: bold;
          line-height: 1;
          cursor: pointer;
          transition: all 0.2s;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
        }

        .controls-collapse-btn:hover {
          background-color: rgba(7, 219, 250, 0.4);
          box-shadow: 0 0 8px rgba(7, 219, 250, 0.4);
        }

        .controls-toggle.collapsed {
          position: absolute;
          top: 10px;
          left: 10px;
          width: 36px;
          height: 36px;
          padding: 0;
          background-color: rgba(10, 20, 30, 0.85);
          backdrop-filter: blur(5px);
          color: #07DBFA;
          border: 1px solid rgba(7, 219, 250, 0.3);
          border-radius: 4px;
          font-size: 20px;
          font-weight: bold;
          cursor: pointer;
          transition: all 0.2s;
          z-index: 10;
          box-shadow: 0 0 15px rgba(7, 219, 250, 0.2);
          display: flex;
          align-items: center;
          justify-content: center;
        }

        .controls-toggle.collapsed:hover {
          background-color: rgba(7, 219, 250, 0.2);
          box-shadow: 0 0 20px rgba(7, 219, 250, 0.4);
        }

        .toggle-icon {
          line-height: 1;
        }

        .control-group {
          display: flex;
          flex-direction: column;
          min-width: 140px;
          margin-bottom: 4px;
        }

        .control-group label {
          margin-bottom: 4px;
          font-size: 11px;
          color: #07DBFA;
          text-transform: uppercase;
          letter-spacing: 0.5px;
        }

        .control-group select,
        .control-group input[type="range"] {
          padding: 4px 8px;
          background-color: rgba(0, 0, 0, 0.4);
          color: #7AFFB2;
          border: 1px solid rgba(7, 219, 250, 0.4);
          border-radius: 3px;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
          font-size: 11px;
          height: 26px;
        }

        .reset-view-button {
          margin-top: 4px;
          background-color: rgba(7, 219, 250, 0.2);
          color: #07DBFA;
          border: 1px solid rgba(7, 219, 250, 0.4);
          padding: 6px 10px;
          border-radius: 3px;
          font-size: 11px;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
          cursor: pointer;
          transition: all 0.2s;
          text-transform: uppercase;
          letter-spacing: 1px;
          width: 100%;
        }

        .reset-view-button:hover {
          background-color: rgba(7, 219, 250, 0.3);
          box-shadow: 0 0 10px rgba(7, 219, 250, 0.4);
        }

        .viewer-container {
          flex-grow: 1;
          height: 100%;
          width: 100%;
          position: relative;
        }
      `}</style>
    </div>
  );
}