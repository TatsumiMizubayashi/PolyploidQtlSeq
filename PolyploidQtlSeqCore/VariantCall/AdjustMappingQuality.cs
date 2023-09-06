namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// Adjust Mapping Quality
    /// </summary>
    internal class AdjustMappingQuality
    {
        /// <summary>
        /// adjust MQの最小値
        /// </summary>
        private const int MINIMUM = 0;

        /// <summary>
        /// adjust MQの最大値
        /// </summary>
        private const int MAXIMUM = 1000;

        /// <summary>
        /// Adjust Mapping Qualityを作成する。
        /// </summary>
        /// <param name="adjustMq">adjustMq</param>
        public AdjustMappingQuality(int adjustMq)
        {
            if (adjustMq < MINIMUM || adjustMq > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(adjustMq));

            Value = adjustMq;
        }

        /// <summary>
        /// Adjust Mapping Qualityを取得する。
        /// </summary>
        internal int Value { get; }
    }
}
