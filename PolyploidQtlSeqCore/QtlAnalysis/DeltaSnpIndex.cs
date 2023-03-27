namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// ΔSNP-index
    /// </summary>
    internal class DeltaSnpIndex
    {
        /// <summary>
        /// ΔSNP-indexを作成する。
        /// </summary>
        /// <param name="value">ΔSNP-index値</param>
        public DeltaSnpIndex(double value)
        {
            Value = value;
        }

        /// <summary>
        /// ΔSNP-index値を取得する。
        /// </summary>
        public double Value { get; }

        /// <summary>
        /// ΔSNP-indexの絶対値を取得する。
        /// </summary>
        /// <returns>ΔSNP-inde絶対値</returns>
        public double Abs()
        {
            return Math.Abs(Value);
        }

        /// <summary>
        /// QTLかどうかを判断する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        /// <returns>QTLならtrue</returns>
        public bool IsQtl(QtlThresholdDeltaSnpIndex threshold)
        {
            return Abs() >= threshold.Value;
        }
    }
}
