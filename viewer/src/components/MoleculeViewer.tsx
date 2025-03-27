import { useRef, useState, useMemo, useCallback } from 'react';
import { Canvas, useFrame, useThree } from '@react-three/fiber';
import { OrbitControls, PerspectiveCamera, Line } from '@react-three/drei';
import * as THREE from 'three';
import { getElementColor, getElementRadius } from '../utils/elementColors';

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
  radius, 
  charge, 
  label, 
  options,
}: { 
  position: [number, number, number]; 
  element: string;
  radius: number;
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
  const displayRadius = useMemo(() => {
    if (options.displayStyle === 'space-filling') {
      return radius * options.atomSizeScale;
    } else if (options.displayStyle === 'ball-and-stick') {
      return radius * 0.5 * options.atomSizeScale;
    } else {
      return radius * 0.3 * options.atomSizeScale;
    }
  }, [radius, options.displayStyle, options.atomSizeScale]);

  // Text label for the atom
  const TextLabel = () => {
    const { camera } = useThree();
    const labelRef = useRef<THREE.Group>(null);
    
    useFrame(() => {
      if (labelRef.current) {
        labelRef.current.lookAt(camera.position);
      }
    });
    
    return options.showLabels ? (
      <group position={position} ref={labelRef}>
        <sprite position={[0, displayRadius + 0.2, 0]} scale={[1, 0.5, 1]}>
          <spriteMaterial attach="material" color="#ffffff">
            <canvasTexture args={[makeTextCanvas(label || element)]} attach="map" />
          </spriteMaterial>
        </sprite>
      </group>
    ) : null;
  };

  // Create a text canvas for the label
  const makeTextCanvas = (text: string): HTMLCanvasElement => {
    const canvas = document.createElement('canvas');
    const ctx = canvas.getContext('2d')!;
    const scale = 5;
    
    canvas.width = text.length * 12 * scale;
    canvas.height = 20 * scale;
    
    ctx.scale(scale, scale);
    ctx.fillStyle = 'rgba(0, 0, 0, 0.6)';
    ctx.fillRect(0, 0, canvas.width / scale, canvas.height / scale);
    
    ctx.font = '12px Arial';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillStyle = '#ffffff';
    ctx.fillText(text, canvas.width / (2 * scale), canvas.height / (2 * scale));
    
    return canvas;
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
    return [
      (start[0] + end[0]) / 2,
      (start[1] + end[1]) / 2,
      (start[2] + end[2]) / 2
    ];
  }, [start, end]);

  // Calculate bond direction and length
  const direction = useMemo(() => {
    const dx = end[0] - start[0];
    const dy = end[1] - start[1];
    const dz = end[2] - start[2];
    const length = Math.sqrt(dx * dx + dy * dy + dz * dz);
    return {
      vector: [dx / length, dy / length, dz / length],
      length,
    };
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
    <mesh
      position={middlePoint as [number, number, number]}
      rotation={[
        Math.atan2(
          Math.sqrt(direction.vector[0] * direction.vector[0] + direction.vector[1] * direction.vector[1]),
          direction.vector[2]
        ),
        0,
        Math.atan2(direction.vector[1], direction.vector[0])
      ]}
    >
      <cylinderGeometry
        args={[bondRadius, bondRadius, direction.length, 12]}
      />
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
          radius={getElementRadius(atom.element)}
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
      <div className="viewer-controls">
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
          gap: 6px 8px;
          padding: 4px 6px;
          background-color: rgba(10, 20, 30, 0.7);
          backdrop-filter: blur(5px);
          color: #7AFFB2;
          z-index: 10;
          border-radius: 4px;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
          border: 1px solid rgba(7, 219, 250, 0.3);
          box-shadow: 0 0 15px rgba(7, 219, 250, 0.2);
          width: 480px;
        }
        
        .control-group {
          display: flex;
          flex-direction: column;
          min-width: 145px;
          margin-bottom: 2px;
        }
        
        .control-group label {
          margin-bottom: 1px;
          font-size: 0.45rem;
          color: #07DBFA;
          text-transform: uppercase;
          letter-spacing: 0.5px;
        }
        
        .control-group select,
        .control-group input[type="range"] {
          padding: 0px 2px;
          background-color: rgba(0, 0, 0, 0.4);
          color: #7AFFB2;
          border: 1px solid rgba(7, 219, 250, 0.4);
          border-radius: 2px;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
          font-size: 0.45rem;
          height: 16px;
        }
        
        .reset-view-button {
          margin-top: 2px;
          background-color: rgba(7, 219, 250, 0.2);
          color: #07DBFA;
          border: 1px solid rgba(7, 219, 250, 0.4);
          padding: 1px 2px;
          border-radius: 2px;
          font-size: 0.45rem;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
          cursor: pointer;
          transition: all 0.2s;
          text-transform: uppercase;
          letter-spacing: 1px;
          width: 100%;
          height: 16px;
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