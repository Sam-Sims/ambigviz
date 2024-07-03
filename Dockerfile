#https://dev.to/rogertorres/first-steps-with-docker-rust-30oi
FROM rust:1.76.0-slim as build
LABEL authors="sammy"

# Update
RUN apt-get update -y && apt-get install -y pkg-config libssl-dev build-essential checkinstall zlib1g-dev

# Create the empty cargo project
RUN USER=root cargo new --bin ambigviz
workdir /ambigviz

# Copy the Cargo.toml and Cargo.lock files
COPY ./Cargo.toml ./Cargo.toml
COPY ./Cargo.lock ./Cargo.lock

# Build the dependencies
RUN cargo build --release
RUN rm src/*.rs
COPY ./src ./src

# Build
RUN rm ./target/release/deps/ambigviz*
RUN cargo build --release

FROM rust:1.76.0-slim
COPY --from=build /ambigviz/target/release/ambigviz /usr/local/bin/ambigviz

ENTRYPOINT ["/usr/local/bin/ambigviz"]