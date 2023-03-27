namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// Bp拡張
    /// </summary>
    internal static class BasepairExtensions
    {
        private const double MBP = 1_000_000.0;

        /// <summary>
        /// bpをMbpに変換する。
        /// </summary>
        /// <param name="bp">bp</param>
        /// <returns>Mbp</returns>
        public static double ToMbp(this int bp)
        {
            return bp / MBP;
        }
    }
}
