# Cerius2 Molecular Viewer

A web-based viewer for Cerius2 molecular model files, focusing on kaolinite mineral surfaces.

## Features

- 3D visualization of molecular structures using Three.js
- Support for Cerius2 Model Structure Interface (.msi) and Cartesian (.car) files
- Interactive rotation, zoom, and pan controls
- Multiple visualization modes (ball-and-stick, space-filling, wireframe)
- Element-based or charge-based coloring
- Detailed atom and bond information
- Structure statistics and measurements

## Getting Started

### Prerequisites

- Node.js (v14 or later)
- npm (v6 or later)

### Installation

1. Clone the repository or download the source code
2. Install dependencies:

```bash
npm install
```

### Development

To start the development server:

```bash
npm run dev
```

This will launch the application at http://localhost:5173.

### Building for Production

To build the application for production:

```bash
npm run build
```

The built files will be in the `dist` directory.

## Usage

1. Upload a Cerius2 model file (.msi, .car, .pdb, or .xyz)
2. The application will parse the file and display the 3D model
3. Use the mouse to interact with the model:
   - Left-click and drag to rotate
   - Right-click and drag to pan
   - Scroll to zoom in/out
4. Use the control panel to adjust display options:
   - Toggle atom labels
   - Change display style
   - Adjust atom and bond sizes
   - Change coloring scheme

## Project Structure

```
/
├── public/             # Static assets
├── src/
│   ├── components/     # React components
│   │   ├── FileUploader.tsx
│   │   ├── MoleculeViewer.tsx
│   │   └── ModelInfo.tsx
│   ├── types/          # TypeScript type definitions
│   │   └── cerius2.ts
│   ├── utils/          # Utility functions
│   │   ├── cerius2Parser.ts
│   │   └── elementColors.ts
│   ├── App.tsx         # Main application component
│   └── main.tsx        # Application entry point
├── index.html          # HTML entry point
└── package.json        # Project configuration
```

## Kaolinite Specifics

This viewer is specially designed for visualizing kaolinite mineral surfaces. Kaolinite has a layered structure with alternating silica tetrahedral and alumina octahedral sheets, with a chemical formula of Al₂Si₂O₅(OH)₄.

The viewer supports:
- Visualizing the layered sheet structure
- Identifying hydroxyl groups on the octahedral sheet
- Examining mineral surface interactions

## Technologies Used

- React
- TypeScript
- Three.js with @react-three/fiber and @react-three/drei
- Vite

## License

This project is licensed under the ISC License.