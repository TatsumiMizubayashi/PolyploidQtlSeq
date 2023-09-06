namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// Mapping Quality最小値
    /// </summary>
    internal class MinmumMappingQuality
    {
        /// <summary>
        /// min-MQの最小値
        /// </summary>
        private const int MINIMUM = 0;

        /// <summary>
        /// min-MQの最大値
        /// </summary>
        private const int MAXIMUM = 60;

        /// <summary>
        /// Mapping Qualityの最小値を作成する。
        /// </summary>
        /// <param name="minMq">MQ最小値</param>
        public MinmumMappingQuality(int minMq)
        {
            if (minMq < MINIMUM || minMq > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(minMq));

            Value = minMq;
        }

        /// <summary>
        /// MQ最小値を取得する。
        /// </summary>
        internal int Value { get; }

    }
}
