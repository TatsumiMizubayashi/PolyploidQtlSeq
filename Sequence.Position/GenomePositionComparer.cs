using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NaturalSort.Extension;

namespace Sequence.Position
{
    /// <summary>
    /// ゲノム位置による比較
    /// </summary>
    public sealed class GenomePositionComparer : IComparer<GenomePosition>
    {
        /// <summary>
        /// シークエンス位置を比較する。
        /// </summary>
        /// <param name="x">位置1</param>
        /// <param name="y">位置2</param>
        /// <returns>比較結果</returns>
        /// <exception cref="NotImplementedException"></exception>
        public int Compare(GenomePosition? x, GenomePosition? y)
        {
            if (x == null) return 0;
            if (y == null) return 0;
            
            // Chr
            var naturalSortComparer = StringComparer.OrdinalIgnoreCase.WithNaturalSort();
            var seqNameComparer = naturalSortComparer.Compare(x.ChrName, y.ChrName);
            if (seqNameComparer != 0) return seqNameComparer;

            // Start
            var startComparer = PositionCompare(x.Start, y.Start);
            if (startComparer != 0) return startComparer;

            // End
            return PositionCompare(x.End, y.End);
        }

        private static int PositionCompare(int x, int y)
        {
            return x == y
                ? 0
                : x > y
                    ? 1
                    : -1;
        }
    }
}
