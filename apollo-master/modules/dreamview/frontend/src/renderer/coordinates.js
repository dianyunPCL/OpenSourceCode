import * as THREE from "three";

export default class Coordinates {
    constructor() {
        this.offset = null;
    }

    isInitialized() {
        return this.offset !== null;
    }

    initialize(x, y) {
        this.offset = {
            x: x,
            y: y
        };
    }

    applyOffset(point, reverse = false) {
        if (this.offset === null) {
            console.error("Offset is not set.");
            return null;
        } else if (isNaN(this.offset.x) || isNaN(this.offset.y)) {
            console.error("Offset contains NaN!");
            return null;
        } else if (isNaN(point.x) || isNaN(point.y)) {
            console.warn("Point contains NaN!");
            return null;
        } else if (!isNaN(point.z)) {
            return new THREE.Vector3(
                reverse ? point.x + this.offset.x : point.x - this.offset.x,
                reverse ? point.y + this.offset.y : point.y - this.offset.y,
                point.z);
        }

        return new THREE.Vector2(reverse ? point.x + this.offset.x : point.x - this.offset.x,
                                 reverse ? point.y + this.offset.y : point.y - this.offset.y);
    }

    applyOffsetToArray(points) {
        return points.map(point => this.applyOffset(point));
    }
}
