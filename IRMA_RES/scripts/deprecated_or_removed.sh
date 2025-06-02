if [ -n "$DOUBLE_LOCAL_PROC" ]; then
    echo "IRMA WARNING: 'DOUBLE_LOCAL_PROC' is set. This variable can no longer be configured manually."
fi

if [ -n "$USE_IRMA_CORE" ]; then
    echo "IRMA WARNING: 'USE_IRMA_CORE' is set. This experimental variable will be removed in a future version."
fi
